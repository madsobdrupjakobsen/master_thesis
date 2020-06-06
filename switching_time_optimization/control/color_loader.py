

colors1 = ["#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352"]

set20_light = '#90e5ca'
set21_light = '#f7b096'

red = [(0.3732256824298347, 0.1963911316160451, 0.1881891580161476),
 (0.5464513648596694, 0.19278226323209027, 0.17637831603229526),
 (0.719677047289504, 0.18917339484813536, 0.16456747404844285),
 (0.8905805459438678, 0.1873484557221581, 0.15433551198257067),
 (0.9151864667435601, 0.2979110598487762, 0.24357298474945543),
 (0.9397923875432526, 0.4084736639753942, 0.33281045751634),
 (0.964398308342945, 0.519036268102012, 0.4220479302832243)]

greens = [(0.20442906574394465, 0.32106113033448674, 0.2376470588235294),
 (0.2089811610918877, 0.4454850698449314, 0.27633986928104576),
 (0.21341022683583233, 0.5665462001794181, 0.31398692810457507),
 (0.2832090221709599, 0.6640369088811995, 0.3801819812892478),
 (0.4146046392413174, 0.7326874279123414, 0.47134691785210825),
 (0.5496501345636293, 0.8032449058054594, 0.5650442137639369)]

colors1 = ["#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352"]
nice_colors = ['#69b3a2']

#plt.style.use('seaborn-dark-palette')

import seaborn as sns
sns.set(style='ticks', palette='Set2')
sns.palplot(sns.color_palette("Paired"))
sns.set(style='ticks', palette="Paired")


col_mod1 = colors1[2] #'C0'
col_mod1_line = '#ccc6a9'
col_mod2 = colors1[4] #'C2'
col_mod2_line = '#91ad64'
col_mod3 = colors1[5] #'C6'
col_idle = colors1[1]

exp_line_style = (0, (5, 1))

paired_col= sns.color_palette("Paired")