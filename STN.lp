\* PSTN Robust Execution LP *\
Maximize
Maximize_time_added_back_to_________constraints_while_decoupling: delta_10_5
 + delta_11_3 + delta_12_15 + delta_13_6 + delta_14_21 + delta_15_12
 + delta_16_4 + delta_17_18 + delta_18_17 + delta_19_9 + delta_1_2
 + delta_20_22 + delta_21_14 + delta_22_20 + delta_2_1 + delta_3_11
 + delta_4_16 + delta_5_10 + delta_6_13 + delta_7_8 + delta_8_7 + delta_9_19
Subject To
_C1: t_0_hi - t_0_lo >= 0
_C10: t_9_hi - t_9_lo >= 0
_C11: t_10_hi - t_10_lo >= 0
_C12: t_11_hi - t_11_lo >= 0
_C13: t_12_hi - t_12_lo >= 0
_C14: t_13_hi - t_13_lo >= 0
_C15: t_14_hi - t_14_lo >= 0
_C16: t_15_hi - t_15_lo >= 0
_C17: t_16_hi - t_16_lo >= 0
_C18: t_17_hi - t_17_lo >= 0
_C19: t_18_hi - t_18_lo >= 0
_C2: t_1_hi - t_1_lo >= 0
_C20: t_19_hi - t_19_lo >= 0
_C21: t_20_hi - t_20_lo >= 0
_C22: t_21_hi - t_21_lo >= 0
_C23: t_22_hi - t_22_lo >= 0
_C24: - t_3_lo + t_4_hi <= 198.24
_C25: t_3_hi - t_4_lo <= -191.69545
_C26: - t_5_lo + t_6_hi <= 23.604
_C27: t_5_hi - t_6_lo <= -8.64
_C28: - t_2_lo + t_9_hi <= 16.761
_C29: t_2_hi - t_9_lo <= -7.89
_C3: t_2_hi - t_2_lo >= 0
_C30: - t_12_lo + t_1_hi <= 15.447
_C31: t_12_hi - t_1_lo <= -8.04
_C32: - t_14_lo + t_7_hi <= 16.25455
_C33: t_14_hi - t_7_lo <= -9.71
_C34: - t_10_lo + t_4_hi <= 577.3344
_C35: t_10_hi - t_4_lo <= 0
_C36: - t_13_lo + t_17_hi <= 12.4925
_C37: t_13_hi - t_17_lo <= -1.82
_C38: - t_19_lo + t_20_hi <= 20.711
_C39: t_19_hi - t_20_lo <= -8.96
_C4: t_3_hi - t_3_lo >= 0
_C40: - t_11_lo + t_21_hi <= 9.02122
_C41: t_11_hi - t_21_lo <= -5.58
_C42: - t_15_lo + t_22_hi <= 203.91185
_C43: t_15_hi - t_22_lo <= -167.08715
_C44: t_16_hi - t_8_lo <= 13.21455
_C45: - t_16_lo + t_8_hi <= -6.67
_C46: - t_10_lo + t_18_hi <= 175.18255
_C47: t_10_hi - t_18_lo <= -141.44395
_C48: t_22_hi - t_3_lo <= 577.3344
_C49: - t_22_lo + t_3_hi <= 0
_C5: t_4_hi - t_4_lo >= 0
_C50: - t_15_lo + t_18_hi <= 577.3344
_C51: t_15_hi - t_18_lo <= 0
_C52: - delta_1_2 - t_1_hi + t_2_hi = 30.0144371
_C53: delta_2_1 - t_1_lo + t_2_lo = 29.1255629
_C54: - delta_7_8 - t_7_hi + t_8_hi = 61.4499162
_C55: delta_8_7 - t_7_lo + t_8_lo = 59.6300838
_C56: - delta_10_5 - t_10_hi + t_5_hi = 50.6296964
_C57: delta_5_10 - t_10_lo + t_5_lo = 49.1303036
_C58: - delta_3_11 + t_11_hi - t_3_hi = 7.94260975
_C59: delta_11_3 + t_11_lo - t_3_lo = 7.70739025
_C6: t_5_hi - t_5_lo >= 0
_C60: - delta_6_13 + t_13_hi - t_6_hi = 36.10969225
_C61: delta_13_6 + t_13_lo - t_6_lo = 35.04030775
_C62: - delta_15_12 + t_12_hi - t_15_hi = 25.0610907
_C63: delta_12_15 + t_12_lo - t_15_lo = 24.3189093
_C64: - delta_16_4 - t_16_hi + t_4_hi = 48.58441095
_C65: delta_4_16 - t_16_lo + t_4_lo = 47.14558905
_C66: - delta_17_18 - t_17_hi + t_18_hi = 50.3251874
_C67: delta_18_17 - t_17_lo + t_18_lo = 48.8348126
_C68: - delta_9_19 + t_19_hi - t_9_hi = 39.7587251
_C69: delta_19_9 + t_19_lo - t_9_lo = 38.5812749
_C7: t_6_hi - t_6_lo >= 0
_C70: - delta_20_22 - t_20_hi + t_22_hi = 53.96406995
_C71: delta_22_20 - t_20_lo + t_22_lo = 52.36593005
_C72: - delta_21_14 + t_14_hi - t_21_hi = 59.63808765
_C73: delta_14_21 + t_14_lo - t_21_lo = 57.87191235
_C8: t_7_hi - t_7_lo >= 0
_C9: t_8_hi - t_8_lo >= 0
Bounds
delta_10_5 <= -2.2460964
delta_11_3 <= -0.35235975
delta_12_15 <= -1.1117907
delta_13_6 <= -1.60194225
delta_14_21 <= -2.64573765
delta_15_12 <= -1.1117907
delta_16_4 <= -2.15536095
delta_17_18 <= -2.2325874
delta_18_17 <= -2.2325874
delta_19_9 <= -1.7638251
delta_1_2 <= -1.3315371
delta_20_22 <= -2.39401995
delta_21_14 <= -2.64573765
delta_22_20 <= -2.39401995
delta_2_1 <= -1.3315371
delta_3_11 <= -0.35235975
delta_4_16 <= -2.15536095
delta_5_10 <= -2.2460964
delta_6_13 <= -1.60194225
delta_7_8 <= -2.7261162
delta_8_7 <= -2.7261162
delta_9_19 <= -1.7638251
t_0_hi = 0
t_0_lo = 0
-inf <= t_10_hi <= inf
-inf <= t_10_lo <= inf
-inf <= t_11_hi <= inf
-inf <= t_11_lo <= inf
-inf <= t_12_hi <= inf
-inf <= t_12_lo <= inf
-inf <= t_13_hi <= inf
-inf <= t_13_lo <= inf
-inf <= t_14_hi <= inf
-inf <= t_14_lo <= inf
-inf <= t_15_hi <= inf
-inf <= t_15_lo <= inf
-inf <= t_16_hi <= inf
-inf <= t_16_lo <= inf
-inf <= t_17_hi <= inf
-inf <= t_17_lo <= inf
-inf <= t_18_hi <= inf
-inf <= t_18_lo <= inf
-inf <= t_19_hi <= inf
-inf <= t_19_lo <= inf
-inf <= t_1_hi <= inf
-inf <= t_1_lo <= inf
-inf <= t_20_hi <= inf
-inf <= t_20_lo <= inf
-inf <= t_21_hi <= inf
-inf <= t_21_lo <= inf
-inf <= t_22_hi <= inf
-inf <= t_22_lo <= inf
-inf <= t_2_hi <= inf
-inf <= t_2_lo <= inf
-inf <= t_3_hi <= inf
-inf <= t_3_lo <= inf
-inf <= t_4_hi <= inf
-inf <= t_4_lo <= inf
-inf <= t_5_hi <= inf
-inf <= t_5_lo <= inf
-inf <= t_6_hi <= inf
-inf <= t_6_lo <= inf
-inf <= t_7_hi <= inf
-inf <= t_7_lo <= inf
-inf <= t_8_hi <= inf
-inf <= t_8_lo <= inf
-inf <= t_9_hi <= inf
-inf <= t_9_lo <= inf
End
