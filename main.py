#%%
#!/usr/bin/env python3
"""
LFAST telescope 
Main integration script for active optics control using multiple wavefront sensing methods.
"""

import sys
import os
from pathlib import Path
import time
import datetime

from arrow import now

# Add git root to path for module imports
GIT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(GIT_ROOT))

from camera_control.camera_control import *
from camera_control.socket_manager import *
from camera_control.high_level_functions import *
from astro_pipeline.onsky_processing import *

socket = SocketClient()


def extended_observation(
    savedir,
    duration_minutes=None,
    num_iterations=None,
    shwfs_pf_interval_sec=30,
    focus_sweep_interval_sec=300,
    focus_sweep_range=(-800, 800),
    focus_sweep_points=9,
    do_tip_tilt_correct=True,
    do_focus_correct=True,
    track_alignment=True,
    plot_pf_output=False,
    pf_exptime=0.01,
    shwfs_exptime=120
):
    """Adaptive focus with periodic SHWFS/PF imaging and focus sweeps."""
    
    setup_socket()
    savedir = normalize_savedir(savedir)
    
    if track_alignment: _alignment_tracking = mirror_alignment() # not sure if it should initialize here or by itself in __main__
    
    # Initialize cameras once
    zwo_cam = ZWOASICamera(ASI_filename='lib\\ASICamera2.dll')
    ids_cam = IDSCamera()
    ids_cam.manual_startup()
    
    base_folder = savedir / f"{datetime.now().strftime('%Y%m%d')}"
    base_folder.mkdir(parents=True, exist_ok=True)
    roi = steer_spot_into_roi(zwo_cam, width=512, height=512, exptime=pf_exptime, nimages=1, savedir=base_folder, roi_factor_safety=0.4)
    pf_exptime = zwo_cam.exptime  # Update exptime so you stop autofocusing
    
    try:
        last_shwfs_pf = last_sweep = time.time()
        start_time = time.time()
        iteration = 0

        
        while True:            
            now = time.time()

            # SHWFS and PF capture with shared timestamp folder
            if now - last_shwfs_pf >= shwfs_pf_interval_sec:
                iteration += 1
                
                # Create single timestamped subfolder for both cameras
                base_subfolder = base_folder / f"{datetime.now().strftime('%H%M%S')}"
                base_subfolder.mkdir(parents=True, exist_ok=True)
                
                # SHWFS capture
                ids_imgs, ids_filenames, _, ids_timestamps = ids_cam.capture_imgs(
                    object_name='shwfs', exptime=shwfs_exptime, nimages=10)
                ids_cam.save_data(save_fits=True, save_bmp=False, imgs=ids_imgs,
                    filenames=ids_filenames, timestamps=ids_timestamps, savedir=base_subfolder)
                print(f"[{iteration}] SHWFS captured")
                
                # Settle time for USB bus to stabilize
                time.sleep(1)
                
                try:
                    extra_header = _alignment_tracking.export_header()
                except NameError:
                    extra_header = None
                              
                # PF capture with retry logic and camera reset
                pf_paths = capture_zwo_with_retry(
                    zwo_cam, object_name='pf', exptime=pf_exptime,
                    nimages=10, savedir=base_subfolder, label=str(iteration), roi=roi,
                    extra_header=extra_header)
                
                if pf_paths:
                    if do_tip_tilt_correct:
                        _correct_tip_tilt(pf_paths[-1])
                    if plot_pf_output:
                        _plot_pf_image(pf_paths[-1])
                
                last_shwfs_pf = time.time()
                now = time.time()
            
            # Focus sweep
            if now - last_sweep >= focus_sweep_interval_sec:
                _focus_sweep_and_correct(zwo_cam, ids_cam, savedir, focus_sweep_range, focus_sweep_points, pf_exptime, nimages=5)
                print(f"[{iteration}] Focus sweep completed")
                last_sweep = time.time()
            
            # Check stopping conditions
            if duration_minutes and time.time() - start_time > duration_minutes * 60:
                break
            if num_iterations and iteration >= num_iterations:
                break
            
            time.sleep(1)
    
    finally:
        ids_cam.manual_shutdown()
        del zwo_cam, ids_cam #investigate this later: seems cranky about zwoasi.__init__.__del__
        socket.close()
        
def steer_spot_into_roi(zwo_cam, cent_x=None, cent_y=None, width=None, height=None, exptime=0.01, nimages=1, roi_factor_safety = 0.8,savedir=None):
    """Steer focus into ROI defined by centroid."""
    current_x, current_y, current_width, current_height = zwo_cam.get_roi()
    current_x, current_y, current_width, current_height = roi_convert_topleft_to_centered(current_x, current_y, current_width, current_height)

    new_x = cent_x if cent_x is not None else current_x
    new_y = cent_y if cent_y is not None else current_y
    new_width = width if width is not None else current_width
    new_height = height if height is not None else current_height
    
    # Define ROI around centroid
    x_min = int(new_x - new_width / 2)
    x_max = int(new_x + new_width / 2)
    y_min = int(new_y - new_height / 2)
    y_max = int(new_y + new_height / 2)
    
    # Capture images and compute centroid offset. 
    temp_savedir = savedir / f"centroiding"
    temp_savedir.mkdir(parents=True, exist_ok=True)
    
    pixel_scale = 0.195  # arcsec/pixel
    gain = 0.5
    iteration = 0
    while True:   
        
        try:
            extra_header = _alignment_tracking.export_header()
        except NameError:
            extra_header = None
        
        pf_paths = capture_zwo_with_retry(zwo_cam, object_name='pf', exptime=exptime, nimages=nimages,
                                        savedir=temp_savedir, extra_header=extra_header)
        exptime = zwo_cam.exptime  # Update exptime so you stop autofocusing

        if pf_paths:
            with fits.open(pf_paths[-1]) as hdul:
                img = hdul[0].data
                #This is assuming offset wrt image center 
                offsets = compute_centroid_offset(img, sigma_threshold=10)

                if offsets is not None:
                    x_offset, y_offset = offsets
                    print(f"ROI iteration {iteration}: Centroid offsets (pixels) = {x_offset:.2f}, {y_offset:.2f}")
                    if np.abs(x_offset) < new_width*roi_factor_safety/2 and np.abs(y_offset) < new_height*roi_factor_safety/2:
                        break
                
                    tip(gain * y_offset * pixel_scale / 2)
                    tilt(gain * x_offset * pixel_scale / 2)
                    time.sleep(0.5)
                iteration += 1

    
    new_x, new_y, new_width, new_height = roi_convert_centered_to_topleft(new_x, new_y, new_width, new_height)
    zwo_cam.set_roi(new_x, new_y, new_width, new_height)
    return [new_x, new_y, new_width, new_height]

    
    
def _focus_sweep_and_correct(zwo_cam, ids_cam, savedir, focus_range, num_points, pf_exptime, nimages=5, exposure_time_lut=None, plot_multiple_sweeps=False):
    """Perform focus sweep, find best focus, apply correction."""
    if exposure_time_lut is  None:
        exposure_time_lut = [-0.01 for x in range(num_points)] 
    backlash = 380  # µm, adjust as needed
    focus_pos = np.linspace(focus_range[0], focus_range[1], num_points)
    focus_deltas = compute_focus_deltas(focus_pos)
    focus_datetime = datetime.now().strftime('%H%M%S')
    
    tracker = {'moved': 0.0}
    try:
        for num, (pos, delta) in enumerate(zip(focus_pos, focus_deltas)):
            focus(delta)
            if num == 0:
                focus(-backlash)
                focus(backlash)
            tracker['moved'] += delta
            time.sleep(0.5)
            subfolder = zwo_cam.create_timestamp_subfolder(savedir,f'FocusSweep_{focus_datetime}')
            
            try:
                extra_header = _alignment_tracking.export_header()
            except NameError:
                extra_header = None

            pf_paths = capture_zwo_with_retry(zwo_cam, object_name='pf', exptime=exposure_time_lut[num],
                                              nimages=nimages, savedir=subfolder, extra_header=extra_header)
            exposure_time_lut[num] = zwo_cam.exptime  # Update LUT with actual exptime used
            _correct_tip_tilt(pf_paths[-1])

    except KeyboardInterrupt:
        print('\nFocus sweep interrupted by user.')
    finally:
        # Always return to original focus position
        if tracker['moved'] != 0:
            print(f'Returning focus to origin (undoing {tracker["moved"]:.1f})...')
            focus(-tracker['moved'])
            
    #Implement Peter's code here for best focus position
    #Need to verify file hierarchy and save these files to a focus_sweep dir.
    subdirs = os.listdir(subfolder.parent)
    reduction = RED(subfolder.parent, subfolder.parent)
    result = reduction.reduce_focus_sweep(focus_arr=focus_pos, subdirs=subdirs, visualize=False, specific_suffix=".zwo.fits", plot_multiple_sweeps=plot_multiple_sweeps)

    if plot_multiple_sweeps:
        return focus_pos, result  # result is fwhm_arr

    optimal_focus = result
    print(f'Optimal Focus: {optimal_focus:.2f} [µm]. Moving there now.')
    focus(-backlash)  # Apply backlash correction before final move
    focus(optimal_focus - tracker['moved'])


def _correct_tip_tilt(fits_path):
    """Correct tip/tilt from centroid."""
    with fits.open(fits_path) as hdul:
        img = hdul[0].data
        offsets = compute_centroid_offset(img, sigma_threshold=10)
        if offsets is None:
            return
        x_offset, y_offset = offsets
        pixel_scale = 0.195  # arcsec/pixel
        gain = 0.5
        tip(gain * y_offset * pixel_scale / 2)
        tilt(gain * x_offset * pixel_scale / 2)
        time.sleep(0.5)

def _plot_pf_image(fits_path):
    """Plot PF image."""
    import matplotlib.pyplot as plt
    with fits.open(fits_path) as hdul:
        plt.imshow(hdul[0].data, cmap='viridis', origin='lower')
        plt.colorbar()
        plt.title('PF')
        plt.show()


def investigate_backlash(
        savedir,
        duration_minutes=30,
        shwfs_pf_interval_sec=30,
        focus_sweep_interval_sec=1,
        focus_sweep_range=(-600, 600),
        focus_sweep_points=9,
        do_tip_tilt_correct=True,
        do_focus_correct=True,
        plot_pf_output=False,
        pf_exptime=-0.1,
        shwfs_exptime=120,
        number_repetitions=5
    ):
    
    setup_socket()
    savedir = normalize_savedir(savedir)
    
    # Initialize cameras once
    zwo_cam = ZWOASICamera(ASI_filename='lib\\ASICamera2.dll')
    ids_cam = IDSCamera()
    ids_cam.manual_startup()
    
    base_folder = savedir / f"{datetime.now().strftime('%Y%m%d')}"
    base_folder.mkdir(parents=True, exist_ok=True)
    roi = steer_spot_into_roi(zwo_cam, width=512, height=512, exptime=pf_exptime, nimages=1, savedir=base_folder, roi_factor_safety=0.4)
    pf_exptime = zwo_cam.exptime  # Update exptime so you stop autofocusing
    
    all_fwhm = []
    for i in range(number_repetitions):
        focus_pos, fwhm_arr = _focus_sweep_and_correct(
            zwo_cam, ids_cam, savedir, focus_sweep_range, focus_sweep_points,
            pf_exptime, nimages=5, plot_multiple_sweeps=True)
        all_fwhm.append(fwhm_arr)
        print(f"[{i}] Focus sweep completed")

    # Overlay all sweeps with sequential colors
    fig, ax = plt.subplots()
    n = len(all_fwhm)
    for i, fwhm_arr in enumerate(all_fwhm):
        color = plt.cm.cool(i / max(n - 1, 1))
        ax.plot(focus_pos, fwhm_arr, 'o-', color=color, label=f'Sweep {i}')
    ax.set_xlabel('Focus Position [µm]')
    ax.set_ylabel('FWHM [pixels]')
    ax.set_title('Backlash Investigation: Repeated Focus Sweeps')
    ax.legend()
    plt.savefig(base_folder / 'backlash_all_sweeps.jpg')
    plt.close()
    
if __name__ == "__main__":
    savedir = '/home/steward/lfast/star_testing/'
    
    #_alignment_tracking = mirror_alignment()
    
    investigate_backlash(
        savedir,
        duration_minutes=30,
        shwfs_pf_interval_sec=30,
        focus_sweep_interval_sec=1,
        focus_sweep_range=(-600, 600),
        focus_sweep_points=9,
        do_tip_tilt_correct=True,
        do_focus_correct=True,
        plot_pf_output=False,
        pf_exptime=-0.1,
        shwfs_exptime=120,
        number_repetitions=5
    )
    
    #_alignment_tracking.reset()    # call if want to reset relative tracking values
    
    extended_observation(
        savedir,
        duration_minutes=30,
        shwfs_pf_interval_sec=30,
        focus_sweep_interval_sec=1,
        focus_sweep_range=(-600, 600),
        focus_sweep_points=9,
        do_tip_tilt_correct=True,
        do_focus_correct=True,
        plot_pf_output=False,
        pf_exptime=-0.1,
        shwfs_exptime=120
    )