

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

blue = [(0.2051672433679354, 0.3085121107266437, 0.38257593233371784),
 (0.2104575163398693, 0.41960784313725485, 0.5694989106753813),
 (0.2174753300012816, 0.5305292836088684, 0.7548225041650647),
 (0.36775599128540315, 0.6244008714596949, 0.8030501089324618),
 (0.5216147635524799, 0.7205074971164935, 0.8524259900038444)]

orange = [(0.9971703191080353, 0.9183391003460207, 0.8395078815840061),
 (0.9943406382160709, 0.8637293348712034, 0.7313802383698578),
 (0.992156862745098, 0.7937254901960784, 0.599769319492503),
 (0.9921568627450981, 0.6933333333333334, 0.4373702422145329),
 (0.992156862745098, 0.5996309111880046, 0.3017916186082278),
 (0.9751787773933104, 0.5020069204152249, 0.17728565936178386),
 (0.937347174163783, 0.40110726643598626, 0.06869665513264145),
 (0.8664821222606689, 0.30366782006920406, 0.015547866205305683),
 (0.7396078431372549, 0.24304498269896194, 0.008289119569396375),
 (0.6083967704728951, 0.19538638985005768, 0.012856593617839307)]

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