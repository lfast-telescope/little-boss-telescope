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

def adaptive_focus_correction(
    savedir,
    duration_minutes=None,
    num_iterations=None,
    shwfs_pf_interval_sec=30,
    focus_sweep_interval_sec=300,
    focus_sweep_range=(-800, 800),
    focus_sweep_points=9,
    do_tip_tilt_correct=True,
    do_focus_correct=True,
    plot_pf_output=False,
    pf_exptime=0.01,
    shwfs_exptime=120
):
    """Adaptive focus with periodic SHWFS/PF imaging and focus sweeps."""
    
    setup_socket()
    if not savedir.endswith('/'):
        savedir += '/'
    savedir = Path(savedir)
    
    # Initialize cameras once
    zwo_cam = ZWOASICamera(ASI_filename='lib\\ASICamera2.dll')
    ids_cam = IDSCamera()
    ids_cam.manual_startup()
    
    base_folder = savedir / f"{datetime.now().strftime('%Y%m%d')}"
    base_folder.mkdir(parents=True, exist_ok=True)
                
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
                
                # CRITICAL: Settle time for USB bus to stabilize
                time.sleep(1)
                              
                # PF capture with retry logic and camera reset
                max_retries = 3
                pf_paths = None
                for attempt in range(max_retries):
                    try:
                        zwo_imgs, zwo_filenames, _, zwo_timestamps = zwo_cam.capture_imgs(
                            object_name='pf', exptime=pf_exptime, nimages=10) 
                        pf_paths = zwo_cam.save_data(save_fits=True, save_bmp=False,
                            imgs=zwo_imgs, filenames=zwo_filenames, timestamps=zwo_timestamps, savedir=base_subfolder)
                        print(f"[{iteration}] PF captured successfully")
                        break
                    except Exception as e:
                        print(f"[{iteration}] PF capture attempt {attempt+1}/{max_retries} failed: {e}")
                        if attempt < max_retries - 1:
                            # Hard reset camera on retry to clear exposure state corruption
                            zwo_cam.camera.stop_video_capture()
                            time.sleep(0.5)
                            print(f"[{iteration}] Camera reset, retrying...")
                            time.sleep(1.5)
                        else:
                            print(f"[{iteration}] PF capture failed after {max_retries} attempts, skipping")
                
                if pf_paths:
                    if do_tip_tilt_correct:
                        _correct_tip_tilt(pf_paths[-1])
                    if plot_pf_output:
                        _plot_pf_image(pf_paths[-1])
                
                last_shwfs_pf = time.time()
            
            # Focus sweep
            if now - last_sweep >= focus_sweep_interval_sec:
                _focus_sweep_and_correct(zwo_cam, ids_cam, savedir, focus_sweep_range, focus_sweep_points)
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
        del zwo_cam, ids_cam
        socket.close()


def _focus_sweep_and_correct(zwo_cam, ids_cam, savedir, focus_range, num_points):
    """Perform focus sweep, find best focus, apply correction."""
    focus_pos = np.linspace(focus_range[0], focus_range[1], num_points)
    focus_deltas = np.insert(np.diff(focus_pos), 0, focus_pos[0])
    
    sharpness = []
    for pos, delta in zip(focus_pos, focus_deltas):
        focus(delta)
        time.sleep(0.5)
        subfolder = zwo_cam.create_timestamp_subfolder(savedir / 'focus_sweep')
        imgs, filenames, _, timestamps = zwo_cam.capture_imgs(
            object_name='sweep', exptime=0.01, nimages=5)
        paths = zwo_cam.save_data(save_fits=True, save_bmp=False,
            imgs=imgs, filenames=filenames, timestamps=timestamps, savedir=subfolder)
            
    #Implement Peter's code here for best focus position


def _correct_tip_tilt(fits_path):
    """Correct tip/tilt from centroid."""
    with fits.open(fits_path) as hdul:
        img = hdul[0].data
        center = np.array(img.shape) / 2
        threshold = np.mean(img) + 10 * np.std(img)
        y_idx, x_idx = np.where(img > threshold)
        if len(x_idx) == 0:
            return
        offset = np.array([np.mean(x_idx) - center[1], np.mean(y_idx) - center[0]])
        offset_arcsec = offset * 0.195
        tip(0.5 * offset_arcsec[1] / 2)
        tilt(0.5 * offset_arcsec[0] / 2)
        time.sleep(0.5)

def _plot_pf_image(fits_path):
    """Plot PF image."""
    import matplotlib.pyplot as plt
    with fits.open(fits_path) as hdul:
        plt.imshow(hdul[0].data, cmap='viridis', origin='lower')
        plt.colorbar()
        plt.title('PF')
        plt.show()

if __name__ == "__main__":
    savedir = '/home/steward/lfast/star_testing/'
    
    adaptive_focus_correction(
        savedir=savedir,
        duration_minutes=30,
        shwfs_pf_interval_sec=1,
        focus_sweep_interval_sec=300,
        focus_sweep_range=(-800, 800),
        focus_sweep_points=9,
        do_tip_tilt_correct=True,
        do_focus_correct=True,
        plot_pf_output=False,
        pf_exptime=0.01,
        shwfs_exptime=120
    )