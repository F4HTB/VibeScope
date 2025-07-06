import sys
import re
import shutil

RENAME_MAP = {
    # AUDIO / ALSA
    "audio_ringbuffer": "audio_ringbuffer",
    "audio_ringbuffer_pos": "audio_ringbuffer_pos",
    "fft_samples": "audio_fft_samples",
    "pcm_handle": "alsa_handle",
    "update_audio_levels": "audio_update_levels",

    # FFT
    "hann_window": "fft_hann_window",
    "FFT_in": "fft_input",
    "FFT_out": "fft_output",
    "fft_plan": "fft_plan",
    "fft_mutex": "fft_mutex",
    "fft_cond": "fft_cond",
    "new_block_ready": "fft_block_ready",

    # SPECTRO
    "col_values": "spectro_band_db",
    "spectro_peak_hold": "spectro_peak_hold",
    "spectro_peak_time": "spectro_peak_time",
    "ref_val": "spectro_ref_val",
    "calibrate_spectro": "spectro_calibrate",
    "band_freqs": "spectro_band_freqs",
    "band_edges": "spectro_band_edges",
    "spectro_tau_peak_hold": "spectro_tau_peak_hold",
    "spectro_tau_peak_fall": "spectro_tau_peak_fall",
    "bar_x": "spectro_bar_x",
    "bar_y": "spectro_bar_y",
    "bar_w": "spectro_bar_w",
    "bar_h": "spectro_bar_h",

    # LSR
    "lr_values": "lsr_rms",
    "lr_peak_values": "lsr_peak",
    "balance_value": "lsr_balance",
    "ref_vu_L": "lsr_ref_vu_L",
    "ref_vu_R": "lsr_ref_vu_R",
    "calibrate_vu": "lsr_calibrate",
    "lr_bar_x": "lsr_bar_x",
    "lr_bar_w": "lsr_bar_w",
    "LR_BAR_Y": "lsr_bar_y",
    "LR_BAR_H": "lsr_bar_h",
    "lr_labels": "lsr_labels",
    "tau_rms": "lsr_tau_rms",
    "tau_peak_release": "lsr_tau_peak_release",

    # UI
    "gradientStops": "ui_gradient_stops",
    "gradTex": "ui_gradient_texture",
    "gradientLUT": "ui_gradient_LUT",
    "circle_cy": "ui_circle_cy",
    "circle_radius": "ui_circle_radius",
    "db_mark": "ui_db_mark",
    "y_mark": "ui_y_mark",
    "NUM_SCALE": "ui_num_scale",

    # Global
    "running": "app_running",
}

def rename_in_file(filename, mapping):
    # Backup
    shutil.copy2(filename, filename + ".bak")

    # Sort by decreasing length (to avoid partial matches)
    sorted_keys = sorted(mapping.keys(), key=len, reverse=True)

    # Pattern: match only variable/function name (not substrings)
    patterns = [re.compile(rf'\b{re.escape(k)}\b') for k in sorted_keys]

    with open(filename, "r") as f:
        code = f.read()

    for pat, k in zip(patterns, sorted_keys):
        v = mapping[k]
        code = pat.sub(v, code)

    with open(filename, "w") as f:
        f.write(code)

    print(f"✅ Renommage appliqué à : {filename} (backup: {filename}.bak)")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage : python rename_vars.py MONFICHIER.c")
        sys.exit(1)

    filename = sys.argv[1]
    rename_in_file(filename, RENAME_MAP)
