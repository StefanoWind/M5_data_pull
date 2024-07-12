# M5_data_pull
Pull data of M5 tower from Y-drive

Run in sequence M5_00_a0>M5_a0_b0>M5_b0_c0.

M5_b0_c0 calculates wind direction both from uncorrected mean velocity components and tilt-corrected ones. Tilt correction is based on [Wilczak et al. 2001](https://link.springer.com/article/10.1023/A:1018966204465) and data from May to August 2022. Wind direction should be minimally affected by tilt anyway. The wind direction are corrected assuming that x is the along-boom component (so 278 deg insted of the canonincal 270 deg), as it appeared from the comparison with lidar data in 2022.
