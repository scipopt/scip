*NAME:         modglob
*ROWS:         291
*COLUMNS:      422
*INTEGER:      98
*NONZERO:      968
*BEST SOLN:    20740508 (opt)
*LP SOLN:      20430947.0
*SOURCE:       Y. Smeers (Univ. of Louvain)
*              Laurence A. Wolsey (Univ. of Louvain)
*              Martin W. P. Savelsbergh (Eindhoven Univ. of Technology)
*APPLICATION:  heating system design, 30 districts, 50 2-way links, 2 products
*COMMENTS:     all integer variables are binary
*              solution reported by Martin W. P. Savelsbergh
*       
NAME          MODGLOB
ROWS
 N  OBJ     
 E  BLDHTXX 
 E  BLGAZXX 
 E  BLDHT15 
 E  BLGAZ15 
 E  DEM01   
 E  DEM02   
 E  DEM03   
 E  DEM04   
 E  DEM05   
 E  DEM06   
 E  DEM07   
 E  DEM08   
 E  DEM09   
 E  DEM10   
 E  DEM11   
 E  DEM12   
 E  DEM13   
 E  DEM14   
 E  DEM16   
 E  DEM17   
 E  DEM18   
 E  DEM19   
 E  DEM20   
 E  DEM21   
 E  DEM22   
 E  DEM23   
 E  DEM24   
 E  DEM25   
 E  DEM26   
 E  DEM27   
 E  DEM28   
 E  DEM29   
 E  DEM30   
 E  DEM31   
 E  DEM15   
 E  CODHT01 
 E  COGAZ01 
 E  CODHT02 
 E  COGAZ02 
 E  CODHT03 
 E  COGAZ03 
 E  CODHT04 
 E  COGAZ04 
 E  CODHT05 
 E  COGAZ05 
 E  CODHT06 
 E  COGAZ06 
 E  CODHT07 
 E  COGAZ07 
 E  CODHT08 
 E  COGAZ08 
 E  CODHT09 
 E  COGAZ09 
 E  CODHT10 
 E  COGAZ10 
 E  CODHT11 
 E  COGAZ11 
 E  CODHT12 
 E  COGAZ12 
 E  CODHT13 
 E  COGAZ13 
 E  CODHT14 
 E  COGAZ14 
 E  CODHT16 
 E  COGAZ16 
 E  CODHT17 
 E  COGAZ17 
 E  CODHT18 
 E  COGAZ18 
 E  CODHT19 
 E  COGAZ19 
 E  CODHT20 
 E  COGAZ20 
 E  CODHT21 
 E  COGAZ21 
 E  CODHT22 
 E  COGAZ22 
 E  CODHT23 
 E  COGAZ23 
 E  CODHT24 
 E  COGAZ24 
 E  CODHT25 
 E  COGAZ25 
 E  CODHT26 
 E  COGAZ26 
 E  CODHT27 
 E  COGAZ27 
 E  CODHT28 
 E  COGAZ28 
 E  CODHT29 
 E  COGAZ29 
 E  CODHT30 
 E  COGAZ30 
 E  CODHT31 
 E  COGAZ31 
 L  LI01DHT 
 L  LI01GAZ 
 L  LI02DHT 
 L  LI02GAZ 
 L  LI03DHT 
 L  LI03GAZ 
 L  LI04DHT 
 L  LI04GAZ 
 L  LI05DHT 
 L  LI05GAZ 
 L  LI06DHT 
 L  LI06GAZ 
 L  LI07DHT 
 L  LI07GAZ 
 L  LI08DHT 
 L  LI08GAZ 
 L  LI09DHT 
 L  LI09GAZ 
 L  LI10DHT 
 L  LI10GAZ 
 L  LI11DHT 
 L  LI11GAZ 
 L  LI12DHT 
 L  LI12GAZ 
 L  LI13DHT 
 L  LI13GAZ 
 L  LI14DHT 
 L  LI14GAZ 
 L  LI15DHT 
 L  LI15GAZ 
 L  LI16DHT 
 L  LI16GAZ 
 L  LI17DHT 
 L  LI17GAZ 
 L  LI18DHT 
 L  LI18GAZ 
 L  LI19DHT 
 L  LI19GAZ 
 L  LI20DHT 
 L  LI20GAZ 
 L  LI21DHT 
 L  LI21GAZ 
 L  LI22DHT 
 L  LI22GAZ 
 L  LI23DHT 
 L  LI23GAZ 
 L  LI24DHT 
 L  LI24GAZ 
 L  LI25DHT 
 L  LI25GAZ 
 L  LI26DHT 
 L  LI26GAZ 
 L  LI27DHT 
 L  LI27GAZ 
 L  LI28DHT 
 L  LI28GAZ 
 L  LI29DHT 
 L  LI29GAZ 
 L  LI30DHT 
 L  LI30GAZ 
 L  LI31DHT 
 L  LI31GAZ 
 L  LI32DHT 
 L  LI32GAZ 
 L  LI33DHT 
 L  LI33GAZ 
 L  LI34DHT 
 L  LI34GAZ 
 L  LI35DHT 
 L  LI35GAZ 
 L  LI36DHT 
 L  LI36GAZ 
 L  LI37DHT 
 L  LI37GAZ 
 L  LI38DHT 
 L  LI38GAZ 
 L  LI39DHT 
 L  LI39GAZ 
 L  LI40DHT 
 L  LI40GAZ 
 L  LI41DHT 
 L  LI41GAZ 
 L  LI42DHT 
 L  LI42GAZ 
 L  LI43DHT 
 L  LI43GAZ 
 L  LI44DHT 
 L  LI44GAZ 
 L  LI45DHT 
 L  LI45GAZ 
 L  LI46DHT 
 L  LI46GAZ 
 L  LI47DHT 
 L  LI47GAZ 
 L  LI48DHT 
 L  LI48GAZ 
 L  LI49DHT 
 L  LI49GAZ 
 L  LJ01DHT 
 L  LJ01GAZ 
 L  LJ02DHT 
 L  LJ02GAZ 
 L  LJ03DHT 
 L  LJ03GAZ 
 L  LJ04DHT 
 L  LJ04GAZ 
 L  LJ05DHT 
 L  LJ05GAZ 
 L  LJ06DHT 
 L  LJ06GAZ 
 L  LJ07DHT 
 L  LJ07GAZ 
 L  LJ08DHT 
 L  LJ08GAZ 
 L  LJ09DHT 
 L  LJ09GAZ 
 L  LJ10DHT 
 L  LJ10GAZ 
 L  LJ11DHT 
 L  LJ11GAZ 
 L  LJ12DHT 
 L  LJ12GAZ 
 L  LJ13DHT 
 L  LJ13GAZ 
 L  LJ14DHT 
 L  LJ14GAZ 
 L  LJ15DHT 
 L  LJ15GAZ 
 L  LJ16DHT 
 L  LJ16GAZ 
 L  LJ17DHT 
 L  LJ17GAZ 
 L  LJ18DHT 
 L  LJ18GAZ 
 L  LJ19DHT 
 L  LJ19GAZ 
 L  LJ20DHT 
 L  LJ20GAZ 
 L  LJ21DHT 
 L  LJ21GAZ 
 L  LJ22DHT 
 L  LJ22GAZ 
 L  LJ23DHT 
 L  LJ23GAZ 
 L  LJ24DHT 
 L  LJ24GAZ 
 L  LJ25DHT 
 L  LJ25GAZ 
 L  LJ26DHT 
 L  LJ26GAZ 
 L  LJ27DHT 
 L  LJ27GAZ 
 L  LJ28DHT 
 L  LJ28GAZ 
 L  LJ29DHT 
 L  LJ29GAZ 
 L  LJ30DHT 
 L  LJ30GAZ 
 L  LJ31DHT 
 L  LJ31GAZ 
 L  LJ32DHT 
 L  LJ32GAZ 
 L  LJ33DHT 
 L  LJ33GAZ 
 L  LJ34DHT 
 L  LJ34GAZ 
 L  LJ35DHT 
 L  LJ35GAZ 
 L  LJ36DHT 
 L  LJ36GAZ 
 L  LJ37DHT 
 L  LJ37GAZ 
 L  LJ38DHT 
 L  LJ38GAZ 
 L  LJ39DHT 
 L  LJ39GAZ 
 L  LJ40DHT 
 L  LJ40GAZ 
 L  LJ41DHT 
 L  LJ41GAZ 
 L  LJ42DHT 
 L  LJ42GAZ 
 L  LJ43DHT 
 L  LJ43GAZ 
 L  LJ44DHT 
 L  LJ44GAZ 
 L  LJ45DHT 
 L  LJ45GAZ 
 L  LJ46DHT 
 L  LJ46GAZ 
 L  LJ47DHT 
 L  LJ47GAZ 
 L  LJ48DHT 
 L  LJ48GAZ 
 L  LJ49DHT 
 L  LJ49GAZ 
COLUMNS
    PRDHTXX   OBJ                432   BLDHTXX           -0.9
    PRDHT15   OBJ                432   BLDHT15           -0.9
    PRGAZXX   OBJ                367   BLGAZXX             -1
    PRGAZ15   OBJ                367   BLGAZ15             -1
    TDHT0102  OBJ            0.65392   CODHT01             -1
    TDHT0102  CODHT02           0.98   LI02DHT              1
    TDHT0104  OBJ           3.224038   CODHT01             -1
    TDHT0104  CODHT04           0.98   LI01DHT              1
    TDHT0201  OBJ            0.65392   CODHT01           0.98
    TDHT0201  CODHT02             -1   LJ02DHT              1
    TDHT0205  OBJ           2.985518   CODHT02             -1
    TDHT0205  CODHT05           0.98   LI03DHT              1
    TDHT0304  OBJ           1.599959   CODHT03             -1
    TDHT0304  CODHT04           0.98   LJ07DHT              1
    TDHT0306  OBJ           4.336906   CODHT03             -1
    TDHT0306  CODHT06           0.98   LI08DHT              1
    TDHT0401  OBJ           3.224038   CODHT01           0.98
    TDHT0401  CODHT04             -1   LJ01DHT              1
    TDHT0403  OBJ           1.599959   CODHT03           0.98
    TDHT0403  CODHT04             -1   LI07DHT              1
    TDHT0405  OBJ           1.135649   CODHT04             -1
    TDHT0405  CODHT05           0.98   LJ04DHT              1
    TDHT0407  OBJ           4.795188   CODHT04             -1
    TDHT0407  CODHT07           0.98   LI06DHT              1
    TDHT0502  OBJ           2.985518   CODHT02           0.98
    TDHT0502  CODHT05             -1   LJ03DHT              1
    TDHT0504  OBJ           1.135649   CODHT04           0.98
    TDHT0504  CODHT05             -1   LI04DHT              1
    TDHT0508  OBJ           5.424316   CODHT05             -1
    TDHT0508  CODHT08           0.98   LI05DHT              1
    TDHT0603  OBJ           4.336906   CODHT03           0.98
    TDHT0603  CODHT06             -1   LJ08DHT              1
    TDHT0607  OBJ           1.938308   CODHT06             -1
    TDHT0607  CODHT07           0.98   LJ09DHT              1
    TDHT0613  OBJ           3.851157   CODHT06             -1
    TDHT0613  CODHT13           0.98   LI15DHT              1
    TDHT0704  OBJ           4.795188   CODHT04           0.98
    TDHT0704  CODHT07             -1   LJ06DHT              1
    TDHT0706  OBJ           1.938308   CODHT06           0.98
    TDHT0706  CODHT07             -1   LI09DHT              1
    TDHT0708  OBJ             1.6951   CODHT07             -1
    TDHT0708  CODHT08           0.98   LI10DHT              1
    TDHT0712  OBJ           4.293357   CODHT07             -1
    TDHT0712  CODHT12           0.98   LI14DHT              1
    TDHT0805  OBJ           5.424316   CODHT05           0.98
    TDHT0805  CODHT08             -1   LJ05DHT              1
    TDHT0807  OBJ             1.6951   CODHT07           0.98
    TDHT0807  CODHT08             -1   LJ10DHT              1
    TDHT0809  OBJ           1.461269   CODHT08             -1
    TDHT0809  CODHT09           0.98   LI11DHT              1
    TDHT0811  OBJ           3.782818   CODHT08             -1
    TDHT0811  CODHT11           0.98   LI13DHT              1
    TDHT0908  OBJ           1.461269   CODHT08           0.98
    TDHT0908  CODHT09             -1   LJ11DHT              1
    TDHT0910  OBJ           3.392878   CODHT09             -1
    TDHT0910  CODHT10           0.98   LI12DHT              1
    TDHT1009  OBJ           3.392878   CODHT09           0.98
    TDHT1009  CODHT10             -1   LJ12DHT              1
    TDHT1011  OBJ           2.129257   CODHT10             -1
    TDHT1011  CODHT11           0.98   LI16DHT              1
    TDHT1015  OBJ           3.449157   CODHT10             -1
    TDHT1015  LI19DHT              1
    TDHT1108  OBJ           3.782818   CODHT08           0.98
    TDHT1108  CODHT11             -1   LJ13DHT              1
    TDHT1110  OBJ           2.129257   CODHT10           0.98
    TDHT1110  CODHT11             -1   LJ16DHT              1
    TDHT1112  OBJ           2.294749   CODHT11             -1
    TDHT1112  CODHT12           0.98   LI17DHT              1
    TDHT1116  OBJ             3.3768   CODHT11             -1
    TDHT1116  CODHT16           0.98   LI20DHT              1
    TDHT1207  OBJ           4.293357   CODHT07           0.98
    TDHT1207  CODHT12             -1   LJ14DHT              1
    TDHT1211  OBJ           2.294749   CODHT11           0.98
    TDHT1211  CODHT12             -1   LJ17DHT              1
    TDHT1213  OBJ           2.068289   CODHT12             -1
    TDHT1213  CODHT13           0.98   LI18DHT              1
    TDHT1217  OBJ           3.107457   CODHT12             -1
    TDHT1217  CODHT17           0.98   LI21DHT              1
    TDHT1306  OBJ           3.851157   CODHT06           0.98
    TDHT1306  CODHT13             -1   LJ15DHT              1
    TDHT1312  OBJ           2.068289   CODHT12           0.98
    TDHT1312  CODHT13             -1   LJ18DHT              1
    TDHT1314  OBJ           1.222749   CODHT13             -1
    TDHT1314  CODHT14           0.98   LI23DHT              1
    TDHT1318  OBJ           2.622378   CODHT13             -1
    TDHT1318  CODHT18           0.98   LI22DHT              1
    TDHT1413  OBJ           1.222749   CODHT13           0.98
    TDHT1413  CODHT14             -1   LJ23DHT              1
    TDHT1419  OBJ           2.252539   CODHT14             -1
    TDHT1419  CODHT19           0.98   LI24DHT              1
    TDHT1611  OBJ             3.3768   CODHT11           0.98
    TDHT1611  CODHT16             -1   LJ20DHT              1
    TDHT1617  OBJ           2.510489   CODHT16             -1
    TDHT1617  CODHT17           0.98   LI28DHT              1
    TDHT1624  OBJ           4.276607   CODHT16             -1
    TDHT1624  CODHT24           0.98   LI27DHT              1
    TDHT1615  OBJ           2.449518   CODHT16             -1
    TDHT1615  LJ26DHT              1
    TDHT1712  OBJ           3.107457   CODHT12           0.98
    TDHT1712  CODHT17             -1   LJ21DHT              1
    TDHT1716  OBJ           2.510489   CODHT16           0.98
    TDHT1716  CODHT17             -1   LJ28DHT              1
    TDHT1718  OBJ           2.576149   CODHT17             -1
    TDHT1718  CODHT18           0.98   LI30DHT              1
    TDHT1723  OBJ           3.762717   CODHT17             -1
    TDHT1723  CODHT23           0.98   LI29DHT              1
    TDHT1813  OBJ           2.622378   CODHT13           0.98
    TDHT1813  CODHT18             -1   LJ22DHT              1
    TDHT1817  OBJ           2.576149   CODHT17           0.98
    TDHT1817  CODHT18             -1   LJ30DHT              1
    TDHT1819  OBJ           2.059578   CODHT18             -1
    TDHT1819  CODHT19           0.98   LI32DHT              1
    TDHT1822  OBJ             2.9078   CODHT18             -1
    TDHT1822  CODHT22           0.98   LI31DHT              1
    TDHT1914  OBJ           2.252539   CODHT14           0.98
    TDHT1914  CODHT19             -1   LJ24DHT              1
    TDHT1918  OBJ           2.059578   CODHT18           0.98
    TDHT1918  CODHT19             -1   LJ32DHT              1
    TDHT1921  OBJ             1.9631   CODHT19             -1
    TDHT1921  CODHT21           0.98   LI33DHT              1
    TDHT2021  OBJ           1.797607   CODHT20             -1
    TDHT2021  CODHT21           0.98   LJ43DHT              1
    TDHT2031  OBJ           1.727258   CODHT20             -1
    TDHT2031  CODHT31           0.98   LI44DHT              1
    TDHT2119  OBJ             1.9631   CODHT19           0.98
    TDHT2119  CODHT21             -1   LJ33DHT              1
    TDHT2120  OBJ           1.797607   CODHT20           0.98
    TDHT2120  CODHT21             -1   LI43DHT              1
    TDHT2122  OBJ           2.842138   CODHT21             -1
    TDHT2122  CODHT22           0.98   LJ41DHT              1
    TDHT2130  OBJ           2.178839   CODHT21             -1
    TDHT2130  CODHT30           0.98   LI42DHT              1
    TDHT2218  OBJ             2.9078   CODHT18           0.98
    TDHT2218  CODHT22             -1   LJ31DHT              1
    TDHT2221  OBJ           2.842138   CODHT21           0.98
    TDHT2221  CODHT22             -1   LI41DHT              1
    TDHT2223  OBJ           4.285317   CODHT22             -1
    TDHT2223  CODHT23           0.98   LJ39DHT              1
    TDHT2229  OBJ           3.573778   CODHT22             -1
    TDHT2229  CODHT29           0.98   LI40DHT              1
    TDHT2317  OBJ           3.762717   CODHT17           0.98
    TDHT2317  CODHT23             -1   LJ29DHT              1
    TDHT2322  OBJ           4.285317   CODHT22           0.98
    TDHT2322  CODHT23             -1   LI39DHT              1
    TDHT2324  OBJ           3.946967   CODHT23             -1
    TDHT2324  CODHT24           0.98   LJ37DHT              1
    TDHT2328  OBJ           4.337578   CODHT23             -1
    TDHT2328  CODHT28           0.98   LI38DHT              1
    TDHT2416  OBJ           4.276607   CODHT16           0.98
    TDHT2416  CODHT24             -1   LJ27DHT              1
    TDHT2423  OBJ           3.946967   CODHT23           0.98
    TDHT2423  CODHT24             -1   LI37DHT              1
    TDHT2425  OBJ           2.714838   CODHT24             -1
    TDHT2425  CODHT25           0.98   LJ35DHT              1
    TDHT2427  OBJ              4.556   CODHT24             -1
    TDHT2427  CODHT27           0.98   LI36DHT              1
    TDHT2524  OBJ           2.714838   CODHT24           0.98
    TDHT2524  CODHT25             -1   LI35DHT              1
    TDHT2526  OBJ           4.397878   CODHT25             -1
    TDHT2526  CODHT26           0.98   LI34DHT              1
    TDHT2515  OBJ           4.060867   CODHT25             -1
    TDHT2515  LJ25DHT              1
    TDHT2625  OBJ           4.397878   CODHT25           0.98
    TDHT2625  CODHT26             -1   LJ34DHT              1
    TDHT2627  OBJ           2.985518   CODHT26             -1
    TDHT2627  CODHT27           0.98   LI45DHT              1
    TDHT2724  OBJ              4.556   CODHT24           0.98
    TDHT2724  CODHT27             -1   LJ36DHT              1
    TDHT2726  OBJ           2.985518   CODHT26           0.98
    TDHT2726  CODHT27             -1   LJ45DHT              1
    TDHT2728  OBJ           5.672216   CODHT27             -1
    TDHT2728  CODHT28           0.98   LI46DHT              1
    TDHT2823  OBJ           4.337578   CODHT23           0.98
    TDHT2823  CODHT28             -1   LJ38DHT              1
    TDHT2827  OBJ           5.672216   CODHT27           0.98
    TDHT2827  CODHT28             -1   LJ46DHT              1
    TDHT2829  OBJ           6.280577   CODHT28             -1
    TDHT2829  CODHT29           0.98   LI47DHT              1
    TDHT2922  OBJ           3.573778   CODHT22           0.98
    TDHT2922  CODHT29             -1   LJ40DHT              1
    TDHT2928  OBJ           6.280577   CODHT28           0.98
    TDHT2928  CODHT29             -1   LJ47DHT              1
    TDHT2930  OBJ           4.728188   CODHT29             -1
    TDHT2930  CODHT30           0.98   LI48DHT              1
    TDHT3021  OBJ           2.178839   CODHT21           0.98
    TDHT3021  CODHT30             -1   LJ42DHT              1
    TDHT3029  OBJ           4.728188   CODHT29           0.98
    TDHT3029  CODHT30             -1   LJ48DHT              1
    TDHT3031  OBJ           2.384528   CODHT30             -1
    TDHT3031  CODHT31           0.98   LI49DHT              1
    TDHT3120  OBJ           1.727258   CODHT20           0.98
    TDHT3120  CODHT31             -1   LJ44DHT              1
    TDHT3130  OBJ           2.384528   CODHT30           0.98
    TDHT3130  CODHT31             -1   LJ49DHT              1
    TDHT1510  OBJ           3.449157   BLDHT15              1
    TDHT1510  CODHT10           0.98   LJ19DHT              1
    TDHT1516  OBJ           2.449518   BLDHT15              1
    TDHT1516  CODHT16           0.98   LI26DHT              1
    TDHT1525  OBJ           4.060867   BLDHT15              1
    TDHT1525  CODHT25           0.98   LI25DHT              1
    TGAZ0102  OBJ            0.35136   COGAZ01             -1
    TGAZ0102  COGAZ02           0.98   LI02GAZ              1
    TGAZ0104  OBJ           1.732319   COGAZ01             -1
    TGAZ0104  COGAZ04           0.98   LI01GAZ              1
    TGAZ0201  OBJ            0.35136   COGAZ01           0.98
    TGAZ0201  COGAZ02             -1   LJ02GAZ              1
    TGAZ0205  OBJ           1.604159   COGAZ02             -1
    TGAZ0205  COGAZ05           0.98   LI03GAZ              1
    TGAZ0304  OBJ            0.85968   COGAZ03             -1
    TGAZ0304  COGAZ04           0.98   LJ07GAZ              1
    TGAZ0306  OBJ           2.330278   COGAZ03             -1
    TGAZ0306  COGAZ06           0.98   LI08GAZ              1
    TGAZ0401  OBJ           1.732319   COGAZ01           0.98
    TGAZ0401  COGAZ04             -1   LJ01GAZ              1
    TGAZ0403  OBJ            0.85968   COGAZ03           0.98
    TGAZ0403  COGAZ04             -1   LI07GAZ              1
    TGAZ0405  OBJ             0.6102   COGAZ04             -1
    TGAZ0405  COGAZ05           0.98   LJ04GAZ              1
    TGAZ0407  OBJ           2.576519   COGAZ04             -1
    TGAZ0407  COGAZ07           0.98   LI06GAZ              1
    TGAZ0502  OBJ           1.604159   COGAZ02           0.98
    TGAZ0502  COGAZ05             -1   LJ03GAZ              1
    TGAZ0504  OBJ             0.6102   COGAZ04           0.98
    TGAZ0504  COGAZ05             -1   LI04GAZ              1
    TGAZ0508  OBJ           2.914559   COGAZ05             -1
    TGAZ0508  COGAZ08           0.98   LI05GAZ              1
    TGAZ0603  OBJ           2.330278   COGAZ03           0.98
    TGAZ0603  COGAZ06             -1   LJ08GAZ              1
    TGAZ0607  OBJ           1.041479   COGAZ06             -1
    TGAZ0607  COGAZ07           0.98   LJ09GAZ              1
    TGAZ0613  OBJ           2.069279   COGAZ06             -1
    TGAZ0613  COGAZ13           0.98   LI15GAZ              1
    TGAZ0704  OBJ           2.576519   COGAZ04           0.98
    TGAZ0704  COGAZ07             -1   LJ06GAZ              1
    TGAZ0706  OBJ           1.041479   COGAZ06           0.98
    TGAZ0706  COGAZ07             -1   LI09GAZ              1
    TGAZ0708  OBJ             0.9108   COGAZ07             -1
    TGAZ0708  COGAZ08           0.98   LI10GAZ              1
    TGAZ0712  OBJ           2.306879   COGAZ07             -1
    TGAZ0712  COGAZ12           0.98   LI14GAZ              1
    TGAZ0805  OBJ           2.914559   COGAZ05           0.98
    TGAZ0805  COGAZ08             -1   LJ05GAZ              1
    TGAZ0807  OBJ             0.9108   COGAZ07           0.98
    TGAZ0807  COGAZ08             -1   LJ10GAZ              1
    TGAZ0809  OBJ            0.78516   COGAZ08             -1
    TGAZ0809  COGAZ09           0.98   LI11GAZ              1
    TGAZ0811  OBJ           2.032559   COGAZ08             -1
    TGAZ0811  COGAZ11           0.98   LI13GAZ              1
    TGAZ0908  OBJ            0.78516   COGAZ08           0.98
    TGAZ0908  COGAZ09             -1   LJ11GAZ              1
    TGAZ0910  OBJ           1.823039   COGAZ09             -1
    TGAZ0910  COGAZ10           0.98   LI12GAZ              1
    TGAZ1009  OBJ           1.823039   COGAZ09           0.98
    TGAZ1009  COGAZ10             -1   LJ12GAZ              1
    TGAZ1011  OBJ           1.144079   COGAZ10             -1
    TGAZ1011  COGAZ11           0.98   LI16GAZ              1
    TGAZ1015  OBJ           1.853279   COGAZ10             -1
    TGAZ1015  LI19GAZ              1
    TGAZ1108  OBJ           2.032559   COGAZ08           0.98
    TGAZ1108  COGAZ11             -1   LJ13GAZ              1
    TGAZ1110  OBJ           1.144079   COGAZ10           0.98
    TGAZ1110  COGAZ11             -1   LJ16GAZ              1
    TGAZ1112  OBJ              1.233   COGAZ11             -1
    TGAZ1112  COGAZ12           0.98   LI17GAZ              1
    TGAZ1116  OBJ             1.8144   COGAZ11             -1
    TGAZ1116  COGAZ16           0.98   LI20GAZ              1
    TGAZ1207  OBJ           2.306879   COGAZ07           0.98
    TGAZ1207  COGAZ12             -1   LJ14GAZ              1
    TGAZ1211  OBJ              1.233   COGAZ11           0.98
    TGAZ1211  COGAZ12             -1   LJ17GAZ              1
    TGAZ1213  OBJ            1.11132   COGAZ12             -1
    TGAZ1213  COGAZ13           0.98   LI18GAZ              1
    TGAZ1217  OBJ           1.669679   COGAZ12             -1
    TGAZ1217  COGAZ17           0.98   LI21GAZ              1
    TGAZ1306  OBJ           2.069279   COGAZ06           0.98
    TGAZ1306  COGAZ13             -1   LJ15GAZ              1
    TGAZ1312  OBJ            1.11132   COGAZ12           0.98
    TGAZ1312  COGAZ13             -1   LJ18GAZ              1
    TGAZ1314  OBJ              0.657   COGAZ13             -1
    TGAZ1314  COGAZ14           0.98   LI23GAZ              1
    TGAZ1318  OBJ           1.409039   COGAZ13             -1
    TGAZ1318  COGAZ18           0.98   LI22GAZ              1
    TGAZ1413  OBJ              0.657   COGAZ13           0.98
    TGAZ1413  COGAZ14             -1   LJ23GAZ              1
    TGAZ1419  OBJ            1.21032   COGAZ14             -1
    TGAZ1419  COGAZ19           0.98   LI24GAZ              1
    TGAZ1611  OBJ             1.8144   COGAZ11           0.98
    TGAZ1611  COGAZ16             -1   LJ20GAZ              1
    TGAZ1617  OBJ           1.348919   COGAZ16             -1
    TGAZ1617  COGAZ17           0.98   LI28GAZ              1
    TGAZ1624  OBJ           2.297879   COGAZ16             -1
    TGAZ1624  COGAZ24           0.98   LI27GAZ              1
    TGAZ1615  OBJ           1.316159   COGAZ16             -1
    TGAZ1615  LJ26GAZ              1
    TGAZ1712  OBJ           1.669679   COGAZ12           0.98
    TGAZ1712  COGAZ17             -1   LJ21GAZ              1
    TGAZ1716  OBJ           1.348919   COGAZ16           0.98
    TGAZ1716  COGAZ17             -1   LJ28GAZ              1
    TGAZ1718  OBJ             1.3842   COGAZ17             -1
    TGAZ1718  COGAZ18           0.98   LI30GAZ              1
    TGAZ1723  OBJ           2.021759   COGAZ17             -1
    TGAZ1723  COGAZ23           0.98   LI29GAZ              1
    TGAZ1813  OBJ           1.409039   COGAZ13           0.98
    TGAZ1813  COGAZ18             -1   LJ22GAZ              1
    TGAZ1817  OBJ             1.3842   COGAZ17           0.98
    TGAZ1817  COGAZ18             -1   LJ30GAZ              1
    TGAZ1819  OBJ           1.106639   COGAZ18             -1
    TGAZ1819  COGAZ19           0.98   LI32GAZ              1
    TGAZ1822  OBJ             1.5624   COGAZ18             -1
    TGAZ1822  COGAZ22           0.98   LI31GAZ              1
    TGAZ1914  OBJ            1.21032   COGAZ14           0.98
    TGAZ1914  COGAZ19             -1   LJ24GAZ              1
    TGAZ1918  OBJ           1.106639   COGAZ18           0.98
    TGAZ1918  COGAZ19             -1   LJ32GAZ              1
    TGAZ1921  OBJ             1.0548   COGAZ19             -1
    TGAZ1921  COGAZ21           0.98   LI33GAZ              1
    TGAZ2021  OBJ           0.965879   COGAZ20             -1
    TGAZ2021  COGAZ21           0.98   LJ43GAZ              1
    TGAZ2031  OBJ           0.928079   COGAZ20             -1
    TGAZ2031  COGAZ31           0.98   LI44GAZ              1
    TGAZ2119  OBJ             1.0548   COGAZ19           0.98
    TGAZ2119  COGAZ21             -1   LJ33GAZ              1
    TGAZ2120  OBJ           0.965879   COGAZ20           0.98
    TGAZ2120  COGAZ21             -1   LI43GAZ              1
    TGAZ2122  OBJ            1.52712   COGAZ21             -1
    TGAZ2122  COGAZ22           0.98   LJ41GAZ              1
    TGAZ2130  OBJ           1.170719   COGAZ21             -1
    TGAZ2130  COGAZ30           0.98   LI42GAZ              1
    TGAZ2218  OBJ             1.5624   COGAZ18           0.98
    TGAZ2218  COGAZ22             -1   LJ31GAZ              1
    TGAZ2221  OBJ            1.52712   COGAZ21           0.98
    TGAZ2221  COGAZ22             -1   LI41GAZ              1
    TGAZ2223  OBJ           2.302559   COGAZ22             -1
    TGAZ2223  COGAZ23           0.98   LJ39GAZ              1
    TGAZ2229  OBJ           1.920239   COGAZ22             -1
    TGAZ2229  COGAZ29           0.98   LI40GAZ              1
    TGAZ2317  OBJ           2.021759   COGAZ17           0.98
    TGAZ2317  COGAZ23             -1   LJ29GAZ              1
    TGAZ2322  OBJ           2.302559   COGAZ22           0.98
    TGAZ2322  COGAZ23             -1   LI39GAZ              1
    TGAZ2324  OBJ           2.120759   COGAZ23             -1
    TGAZ2324  COGAZ24           0.98   LJ37GAZ              1
    TGAZ2328  OBJ           2.330639   COGAZ23             -1
    TGAZ2328  COGAZ28           0.98   LI38GAZ              1
    TGAZ2416  OBJ           2.297879   COGAZ16           0.98
    TGAZ2416  COGAZ24             -1   LJ27GAZ              1
    TGAZ2423  OBJ           2.120759   COGAZ23           0.98
    TGAZ2423  COGAZ24             -1   LI37GAZ              1
    TGAZ2425  OBJ           1.458719   COGAZ24             -1
    TGAZ2425  COGAZ25           0.98   LJ35GAZ              1
    TGAZ2427  OBJ              2.448   COGAZ24             -1
    TGAZ2427  COGAZ27           0.98   LI36GAZ              1
    TGAZ2524  OBJ           1.458719   COGAZ24           0.98
    TGAZ2524  COGAZ25             -1   LI35GAZ              1
    TGAZ2526  OBJ           2.363039   COGAZ25             -1
    TGAZ2526  COGAZ26           0.98   LI34GAZ              1
    TGAZ2515  OBJ           2.181959   COGAZ25             -1
    TGAZ2515  LJ25GAZ              1
    TGAZ2625  OBJ           2.363039   COGAZ25           0.98
    TGAZ2625  COGAZ26             -1   LJ34GAZ              1
    TGAZ2627  OBJ           1.604159   COGAZ26             -1
    TGAZ2627  COGAZ27           0.98   LI45GAZ              1
    TGAZ2724  OBJ              2.448   COGAZ24           0.98
    TGAZ2724  COGAZ27             -1   LJ36GAZ              1
    TGAZ2726  OBJ           1.604159   COGAZ26           0.98
    TGAZ2726  COGAZ27             -1   LJ45GAZ              1
    TGAZ2728  OBJ           3.047759   COGAZ27             -1
    TGAZ2728  COGAZ28           0.98   LI46GAZ              1
    TGAZ2823  OBJ           2.330639   COGAZ23           0.98
    TGAZ2823  COGAZ28             -1   LJ38GAZ              1
    TGAZ2827  OBJ           3.047759   COGAZ27           0.98
    TGAZ2827  COGAZ28             -1   LJ46GAZ              1
    TGAZ2829  OBJ            3.37464   COGAZ28             -1
    TGAZ2829  COGAZ29           0.98   LI47GAZ              1
    TGAZ2922  OBJ           1.920239   COGAZ22           0.98
    TGAZ2922  COGAZ29             -1   LJ40GAZ              1
    TGAZ2928  OBJ            3.37464   COGAZ28           0.98
    TGAZ2928  COGAZ29             -1   LJ47GAZ              1
    TGAZ2930  OBJ            2.54052   COGAZ29             -1
    TGAZ2930  COGAZ30           0.98   LI48GAZ              1
    TGAZ3021  OBJ           1.170719   COGAZ21           0.98
    TGAZ3021  COGAZ30             -1   LJ42GAZ              1
    TGAZ3029  OBJ            2.54052   COGAZ29           0.98
    TGAZ3029  COGAZ30             -1   LJ48GAZ              1
    TGAZ3031  OBJ            1.28124   COGAZ30             -1
    TGAZ3031  COGAZ31           0.98   LI49GAZ              1
    TGAZ3120  OBJ           0.928079   COGAZ20           0.98
    TGAZ3120  COGAZ31             -1   LJ44GAZ              1
    TGAZ3130  OBJ            1.28124   COGAZ30           0.98
    TGAZ3130  COGAZ31             -1   LJ49GAZ              1
    TGAZ1510  OBJ           1.853279   BLGAZ15              1
    TGAZ1510  COGAZ10           0.98   LJ19GAZ              1
    TGAZ1516  OBJ           1.316159   BLGAZ15              1
    TGAZ1516  COGAZ16           0.98   LI26GAZ              1
    TGAZ1525  OBJ           2.181959   BLGAZ15              1
    TGAZ1525  COGAZ25           0.98   LI25GAZ              1
    MARK0000  'MARKER'                 'INTORG'
    D01DHT    OBJ        15061.55078   LI01DHT         -17360
    D01DHT    LJ01DHT         -17360
    D01GAZ    OBJ        9479.636719   LI01GAZ         -17360
    D01GAZ    LJ01GAZ         -17360
    D02DHT    OBJ        3054.878418   LI02DHT         -17360
    D02DHT    LJ02DHT         -17360
    D02GAZ    OBJ        1922.719482   LI02GAZ         -17360
    D02GAZ    LJ02GAZ         -17360
    D03DHT    OBJ        13947.26953   LI03DHT         -17360
    D03DHT    LJ03DHT         -17360
    D03GAZ    OBJ          8778.3125   LI03GAZ         -17360
    D03GAZ    LJ03GAZ         -17360
    D04DHT    OBJ        5305.347656   LI04DHT         -17360
    D04DHT    LJ04DHT         -17360
    D04GAZ    OBJ        3339.149414   LI04GAZ         -17360
    D04GAZ    LJ04GAZ         -17360
    D05DHT    OBJ        25340.46484   LI05DHT         -17360
    D05DHT    LJ05DHT         -17360
    D05GAZ    OBJ        15949.11328   LI05GAZ         -17360
    D05GAZ    LJ05GAZ         -17360
    D06DHT    OBJ        22401.39844   LI06DHT         -17360
    D06DHT    LJ06DHT         -17360
    D06GAZ    OBJ        14099.28516   LI06GAZ         -17360
    D06GAZ    LJ06GAZ         -17360
    D07DHT    OBJ        7474.433594   LI07DHT         -17360
    D07DHT    LJ07DHT         -17360
    D07GAZ    OBJ        4704.355469   LI07GAZ         -17360
    D07GAZ    LJ07GAZ         -17360
    D08DHT    OBJ        20260.47266   LI08DHT         -17360
    D08DHT    LJ08DHT         -17360
    D08GAZ    OBJ        12751.80078   LI08GAZ         -17360
    D08GAZ    LJ08GAZ         -17360
    D09DHT    OBJ        9055.078125   LI09DHT         -17360
    D09DHT    LJ09DHT         -17360
    D09GAZ    OBJ        5699.203125   LI09GAZ         -17360
    D09GAZ    LJ09GAZ         -17360
    D10DHT    OBJ        7918.894531   LI10DHT         -17360
    D10DHT    LJ10DHT         -17360
    D10GAZ    OBJ        4984.097656   LI10GAZ         -17360
    D10GAZ    LJ10GAZ         -17360
    D11DHT    OBJ        6826.523438   LI11DHT         -17360
    D11DHT    LJ11DHT         -17360
    D11GAZ    OBJ        4296.566406   LI11GAZ         -17360
    D11GAZ    LJ11GAZ         -17360
    D12DHT    OBJ        15850.30859   LI12DHT         -17360
    D12DHT    LJ12DHT         -17360
    D12GAZ    OBJ        9976.074219   LI12GAZ         -17360
    D12GAZ    LJ12GAZ         -17360
    D13DHT    OBJ        17671.96484   LI13DHT         -17360
    D13DHT    LJ13DHT         -17360
    D13GAZ    OBJ        11122.61328   LI13GAZ         -17360
    D13GAZ    LJ13GAZ         -17360
    D14DHT    OBJ        20057.02344   LI14DHT         -17360
    D14DHT    LJ14DHT         -17360
    D14GAZ    OBJ        12623.75391   LI14GAZ         -17360
    D14GAZ    LJ14GAZ         -17360
    D15DHT    OBJ        17991.22656   LI15DHT         -17360
    D15DHT    LJ15DHT         -17360
    D15GAZ    OBJ        11323.55078   LI15GAZ         -17360
    D15GAZ    LJ15GAZ         -17360
    D16DHT    OBJ        9947.128906   LI16DHT         -17360
    D16DHT    LJ16DHT         -17360
    D16GAZ    OBJ        6260.652344   LI16GAZ         -17360
    D16GAZ    LJ16GAZ         -17360
    D17DHT    OBJ        10720.24219   LI17DHT         -17360
    D17DHT    LJ17DHT         -17360
    D17GAZ    OBJ        6747.246094   LI17GAZ         -17360
    D17GAZ    LJ17GAZ         -17360
    D18DHT    OBJ        9662.304688   LI18DHT         -17360
    D18DHT    LJ18DHT         -17360
    D18GAZ    OBJ        6081.386719   LI18GAZ         -17360
    D18GAZ    LJ18GAZ         -17360
    D19DHT    OBJ        16113.22656   LI19DHT         -17360
    D19DHT    LJ19DHT         -17360
    D19GAZ    OBJ        10141.55078   LI19GAZ         -17360
    D19GAZ    LJ19GAZ         -17360
    D20DHT    OBJ        15775.19141   LI20DHT         -17360
    D20DHT    LJ20DHT         -17360
    D20GAZ    OBJ        9928.796875   LI20GAZ         -17360
    D20GAZ    LJ20GAZ         -17360
    D21DHT    OBJ        14516.92578   LI21DHT         -17360
    D21DHT    LJ21DHT         -17360
    D21GAZ    OBJ        9136.851563   LI21GAZ         -17360
    D21GAZ    LJ21GAZ         -17360
    D22DHT    OBJ        12250.80859   LI22DHT         -17360
    D22DHT    LJ22DHT         -17360
    D22GAZ    OBJ        7710.574219   LI22GAZ         -17360
    D22GAZ    LJ22GAZ         -17360
    D23DHT    OBJ        5712.246094   LI23DHT         -17360
    D23DHT    LJ23DHT         -17360
    D23GAZ    OBJ        3595.249268   LI23GAZ         -17360
    D23GAZ    LJ23GAZ         -17360
    D24DHT    OBJ        10523.05078   LI24DHT         -17360
    D24DHT    LJ24DHT         -17360
    D24GAZ    OBJ        6623.136719   LI24GAZ         -17360
    D24GAZ    LJ24GAZ         -17360
    D25DHT    OBJ        18970.91797   LI25DHT         -17360
    D25DHT    LJ25DHT         -17360
    D25GAZ    OBJ        11940.16406   LI25GAZ         -17360
    D25GAZ    LJ25GAZ         -17360
    D26DHT    OBJ        11443.26953   LI26DHT         -17360
    D26DHT    LJ26DHT         -17360
    D26GAZ    OBJ          7202.3125   LI26GAZ         -17360
    D26GAZ    LJ26GAZ         -17360
    D27DHT    OBJ        19978.77344   LI27DHT         -17360
    D27DHT    LJ27DHT         -17360
    D27GAZ    OBJ        12574.50391   LI27GAZ         -17360
    D27GAZ    LJ27GAZ         -17360
    D28DHT    OBJ        11728.10156   LI28DHT         -17360
    D28DHT    LJ28DHT         -17360
    D28GAZ    OBJ        7381.585938   LI28GAZ         -17360
    D28GAZ    LJ28GAZ         -17360
    D29DHT    OBJ        17578.06641   LI29DHT         -17360
    D29DHT    LJ29DHT         -17360
    D29GAZ    OBJ        11063.51172   LI29GAZ         -17360
    D29GAZ    LJ29GAZ         -17360
    D30DHT    OBJ        12034.84375   LI30DHT         -17360
    D30DHT    LJ30DHT         -17360
    D30GAZ    OBJ        7574.648438   LI30GAZ         -17360
    D30GAZ    LJ30GAZ         -17360
    D31DHT    OBJ        13584.19141   LI31DHT         -17360
    D31DHT    LJ31DHT         -17360
    D31GAZ    OBJ        8549.796875   LI31GAZ         -17360
    D31GAZ    LJ31GAZ         -17360
    D32DHT    OBJ        9621.609375   LI32DHT         -17360
    D32DHT    LJ32DHT         -17360
    D32GAZ    OBJ        6055.773438   LI32GAZ         -17360
    D32GAZ    LJ32GAZ         -17360
    D33DHT    OBJ        9170.894531   LI33DHT         -17360
    D33DHT    LJ33DHT         -17360
    D33GAZ    OBJ        5772.097656   LI33GAZ         -17360
    D33GAZ    LJ33GAZ         -17360
    D34DHT    OBJ        20545.30859   LI34DHT         -17360
    D34DHT    LJ34DHT         -17360
    D34GAZ    OBJ        12931.07422   LI34GAZ         -17360
    D34GAZ    LJ34GAZ         -17360
    D35DHT    OBJ           12682.75   LI35DHT         -17360
    D35DHT    LJ35DHT         -17360
    D35GAZ    OBJ          7982.4375   LI35GAZ         -17360
    D35GAZ    LJ35GAZ         -17360
    D36DHT    OBJ        21283.98828   LI36DHT         -17360
    D36DHT    LJ36DHT         -17360
    D36GAZ    OBJ        13395.99609   LI36GAZ         -17360
    D36GAZ    LJ36GAZ         -17360
    D37DHT    OBJ        18438.81641   LI37DHT         -17360
    D37DHT    LJ37DHT         -17360
    D37GAZ    OBJ        11605.26172   LI37GAZ         -17360
    D37GAZ    LJ37GAZ         -17360
    D38DHT    OBJ        20263.60547   LI38DHT         -17360
    D38DHT    LJ38DHT         -17360
    D38GAZ    OBJ        12753.77344   LI38GAZ         -17360
    D38GAZ    LJ38GAZ         -17360
    D39DHT    OBJ        20019.46484   LI39DHT         -17360
    D39DHT    LJ39DHT         -17360
    D39GAZ    OBJ        12600.11328   LI39GAZ         -17360
    D39GAZ    LJ39GAZ         -17360
    D40DHT    OBJ        16695.41016   LI40DHT         -17360
    D40DHT    LJ40DHT         -17360
    D40GAZ    OBJ        10507.97266   LI40GAZ         -17360
    D40GAZ    LJ40GAZ         -17360
    D41DHT    OBJ        13277.45313   LI41DHT         -17360
    D41DHT    LJ41DHT         -17360
    D41GAZ    OBJ        8356.734375   LI41GAZ         -17360
    D41GAZ    LJ41GAZ         -17360
    D42DHT    OBJ        10178.75391   LI42DHT         -17360
    D42DHT    LJ42DHT         -17360
    D42GAZ    OBJ          6406.4375   LI42GAZ         -17360
    D42GAZ    LJ42GAZ         -17360
    D43DHT    OBJ        8397.777344   LI43DHT         -17360
    D43DHT    LJ43DHT         -17360
    D43GAZ    OBJ        5285.503906   LI43GAZ         -17360
    D43GAZ    LJ43GAZ         -17360
    D44DHT    OBJ        8069.128906   LI44DHT         -17360
    D44DHT    LJ44DHT         -17360
    D44GAZ    OBJ        5078.652344   LI44GAZ         -17360
    D44GAZ    LJ44GAZ         -17360
    D45DHT    OBJ        13947.26953   LI45DHT         -17360
    D45DHT    LJ45DHT         -17360
    D45GAZ    OBJ          8778.3125   LI45GAZ         -17360
    D45GAZ    LJ45GAZ         -17360
    D46DHT    OBJ         26498.5625   LI46DHT         -17360
    D46DHT    LJ46DHT         -17360
    D46GAZ    OBJ        16678.01172   LI46GAZ         -17360
    D46GAZ    LJ46GAZ         -17360
    D47DHT    OBJ        29340.60547   LI47DHT         -17360
    D47DHT    LJ47DHT         -17360
    D47GAZ    OBJ        18466.77344   LI47GAZ         -17360
    D47GAZ    LJ47GAZ         -17360
    D48DHT    OBJ        22088.39844   LI48DHT         -17360
    D48DHT    LJ48DHT         -17360
    D48GAZ    OBJ        13902.28516   LI48GAZ         -17360
    D48GAZ    LJ48GAZ         -17360
    D49DHT    OBJ        11139.66016   LI49DHT         -17360
    D49DHT    LJ49DHT         -17360
    D49GAZ    OBJ        7011.226563   LI49GAZ         -17360
    D49GAZ    LJ49GAZ         -17360
    MARK0001  'MARKER'                 'INTEND'
    CDHT01    OBJ         935.779785   DEM01                1
    CDHT01    CODHT01      -1.197604
    CDHT02    OBJ          590.19458   DEM02                1
    CDHT02    CODHT02      -1.197604
    CDHT03    OBJ         586.431152   DEM03                1
    CDHT03    CODHT03      -1.197604
    CDHT04    OBJ         586.384521   DEM04                1
    CDHT04    CODHT04      -1.197604
    CDHT05    OBJ         585.222412   DEM05                1
    CDHT05    CODHT05      -1.197604
    CDHT06    OBJ           584.8667   DEM06                1
    CDHT06    CODHT06      -1.197604
    CDHT07    OBJ         584.592285   DEM07                1
    CDHT07    CODHT07      -1.197604
    CDHT08    OBJ         586.196533   DEM08                1
    CDHT08    CODHT08      -1.197604
    CDHT09    OBJ         607.073486   DEM09                1
    CDHT09    CODHT09      -1.197604
    CDHT10    OBJ         593.749756   DEM10                1
    CDHT10    CODHT10      -1.197604
    CDHT11    OBJ          585.37793   DEM11                1
    CDHT11    CODHT11      -1.197604
    CDHT12    OBJ         584.541016   DEM12                1
    CDHT12    CODHT12      -1.197604
    CDHT13    OBJ         585.455078   DEM13                1
    CDHT13    CODHT13      -1.197604
    CDHT14    OBJ         587.869385   DEM14                1
    CDHT14    CODHT14      -1.197604
    CDHT16    OBJ         584.969482   DEM16                1
    CDHT16    CODHT16      -1.197604
    CDHT17    OBJ         584.601318   DEM17                1
    CDHT17    CODHT17      -1.197604
    CDHT18    OBJ         585.421143   DEM18                1
    CDHT18    CODHT18      -1.197604
    CDHT19    OBJ         586.946777   DEM19                1
    CDHT19    CODHT19      -1.197604
    CDHT20    OBJ         592.529541   DEM20                1
    CDHT20    CODHT20      -1.197604
    CDHT21    OBJ         585.983643   DEM21                1
    CDHT21    CODHT21      -1.197604
    CDHT22    OBJ         585.921387   DEM22                1
    CDHT22    CODHT22      -1.197604
    CDHT23    OBJ         586.138184   DEM23                1
    CDHT23    CODHT23      -1.197604
    CDHT24    OBJ         585.365967   DEM24                1
    CDHT24    CODHT24      -1.197604
    CDHT25    OBJ           593.9563   DEM25                1
    CDHT25    CODHT25      -1.197604
    CDHT26    OBJ         586.775391   DEM26                1
    CDHT26    CODHT26      -1.197604
    CDHT27    OBJ         590.526367   DEM27                1
    CDHT27    CODHT27      -1.197604
    CDHT28    OBJ         588.990479   DEM28                1
    CDHT28    CODHT28      -1.197604
    CDHT29    OBJ         628.276367   DEM29                1
    CDHT29    CODHT29      -1.197604
    CDHT30    OBJ         589.911133   DEM30                1
    CDHT30    CODHT30      -1.197604
    CDHT31    OBJ         590.110107   DEM31                1
    CDHT31    CODHT31      -1.197604
    CDHT15    OBJ         585.704346   BLDHT15       1.197604
    CDHT15    DEM15                1
    CGAZ01    OBJ         994.519775   DEM01                1
    CGAZ01    COGAZ01      -1.503759
    CGAZ02    OBJ         764.129639   DEM02                1
    CGAZ02    COGAZ02      -1.503759
    CGAZ03    OBJ          761.62085   DEM03                1
    CGAZ03    COGAZ03      -1.503759
    CGAZ04    OBJ           761.5896   DEM04                1
    CGAZ04    COGAZ04      -1.503759
    CGAZ05    OBJ         760.814697   DEM05                1
    CGAZ05    COGAZ05      -1.503759
    CGAZ06    OBJ         760.577881   DEM06                1
    CGAZ06    COGAZ06      -1.503759
    CGAZ07    OBJ         760.394287   DEM07                1
    CGAZ07    COGAZ07      -1.503759
    CGAZ08    OBJ         761.463867   DEM08                1
    CGAZ08    COGAZ08      -1.503759
    CGAZ09    OBJ          775.38208   DEM09                1
    CGAZ09    COGAZ09      -1.503759
    CGAZ10    OBJ         766.499756   DEM10                1
    CGAZ10    COGAZ10      -1.503759
    CGAZ11    OBJ         760.918457   DEM11                1
    CGAZ11    COGAZ11      -1.503759
    CGAZ12    OBJ         760.360107   DEM12                1
    CGAZ12    COGAZ12      -1.503759
    CGAZ13    OBJ         760.969971   DEM13                1
    CGAZ13    COGAZ13      -1.503759
    CGAZ14    OBJ         762.579346   DEM14                1
    CGAZ14    COGAZ14      -1.503759
    CGAZ16    OBJ          760.64624   DEM16                1
    CGAZ16    COGAZ16      -1.503759
    CGAZ17    OBJ              760.4   DEM17                1
    CGAZ17    COGAZ17      -1.503759
    CGAZ18    OBJ         760.947021   DEM18                1
    CGAZ18    COGAZ18      -1.503759
    CGAZ19    OBJ         761.964355   DEM19                1
    CGAZ19    COGAZ19      -1.503759
    CGAZ20    OBJ         765.686279   DEM20                1
    CGAZ20    COGAZ20      -1.503759
    CGAZ21    OBJ         761.322021   DEM21                1
    CGAZ21    COGAZ21      -1.503759
    CGAZ22    OBJ         761.280518   DEM22                1
    CGAZ22    COGAZ22      -1.503759
    CGAZ23    OBJ         761.425293   DEM23                1
    CGAZ23    COGAZ23      -1.503759
    CGAZ24    OBJ           760.9104   DEM24                1
    CGAZ24    COGAZ24      -1.503759
    CGAZ25    OBJ         766.636475   DEM25                1
    CGAZ25    COGAZ25      -1.503759
    CGAZ26    OBJ         761.849609   DEM26                1
    CGAZ26    COGAZ26      -1.503759
    CGAZ27    OBJ          764.35083   DEM27                1
    CGAZ27    COGAZ27      -1.503759
    CGAZ28    OBJ          763.32666   DEM28                1
    CGAZ28    COGAZ28      -1.503759
    CGAZ29    OBJ         789.517578   DEM29                1
    CGAZ29    COGAZ29      -1.503759
    CGAZ30    OBJ         763.940186   DEM30                1
    CGAZ30    COGAZ30      -1.503759
    CGAZ31    OBJ         764.072754   DEM31                1
    CGAZ31    COGAZ31      -1.503759
    CGAZ15    OBJ         761.135986   BLGAZ15       1.503759
    CGAZ15    DEM15                1
    CEHT01    OBJ              40996   DEM01                1
    CEHT02    OBJ              40996   DEM02                1
    CEHT03    OBJ        40995.97266   DEM03                1
    CEHT04    OBJ        40995.97656   DEM04                1
    CEHT05    OBJ        40995.99609   DEM05                1
    CEHT06    OBJ        40995.99219   DEM06                1
    CEHT07    OBJ        40995.99609   DEM07                1
    CEHT08    OBJ        40995.99219   DEM08                1
    CEHT09    OBJ        40995.99609   DEM09                1
    CEHT10    OBJ        40995.99609   DEM10                1
    CEHT11    OBJ        40995.97266   DEM11                1
    CEHT12    OBJ        40995.99609   DEM12                1
    CEHT13    OBJ        40995.98438   DEM13                1
    CEHT14    OBJ        40995.99609   DEM14                1
    CEHT16    OBJ        40995.98828   DEM16                1
    CEHT17    OBJ        40995.99609   DEM17                1
    CEHT18    OBJ        40995.99219   DEM18                1
    CEHT19    OBJ        40995.99609   DEM19                1
    CEHT20    OBJ              40996   DEM20                1
    CEHT21    OBJ        40995.99609   DEM21                1
    CEHT22    OBJ        40995.97656   DEM22                1
    CEHT23    OBJ        40995.98438   DEM23                1
    CEHT24    OBJ        40995.98047   DEM24                1
    CEHT25    OBJ        40995.98828   DEM25                1
    CEHT26    OBJ        40995.98828   DEM26                1
    CEHT27    OBJ        40995.99219   DEM27                1
    CEHT28    OBJ        40995.98828   DEM28                1
    CEHT29    OBJ              40996   DEM29                1
    CEHT30    OBJ        40995.99609   DEM30                1
    CEHT31    OBJ        40995.99609   DEM31                1
    CEHT15    OBJ        40995.98047   DEM15                1
    CNEH01    OBJ           2184.375   DEM01                1
    CNEH02    OBJ           2184.375   DEM02                1
    CNEH03    OBJ           2184.375   DEM03                1
    CNEH04    OBJ        2184.374268   DEM04                1
    CNEH05    OBJ        2184.374512   DEM05                1
    CNEH06    OBJ           2184.375   DEM06                1
    CNEH07    OBJ        2184.374512   DEM07                1
    CNEH08    OBJ        2184.374512   DEM08                1
    CNEH09    OBJ        2184.374512   DEM09                1
    CNEH10    OBJ        2184.373779   DEM10                1
    CNEH11    OBJ        2184.374268   DEM11                1
    CNEH12    OBJ        2184.374268   DEM12                1
    CNEH13    OBJ        2184.374268   DEM13                1
    CNEH14    OBJ        2184.374512   DEM14                1
    CNEH16    OBJ        2184.374512   DEM16                1
    CNEH17    OBJ        2184.374512   DEM17                1
    CNEH18    OBJ        2184.374268   DEM18                1
    CNEH19    OBJ        2184.374512   DEM19                1
    CNEH20    OBJ           2184.375   DEM20                1
    CNEH21    OBJ        2184.374268   DEM21                1
    CNEH22    OBJ        2184.374512   DEM22                1
    CNEH23    OBJ        2184.374512   DEM23                1
    CNEH24    OBJ        2184.374512   DEM24                1
    CNEH25    OBJ        2184.374268   DEM25                1
    CNEH26    OBJ        2184.374512   DEM26                1
    CNEH27    OBJ        2184.374512   DEM27                1
    CNEH28    OBJ        2184.374512   DEM28                1
    CNEH29    OBJ           2184.375   DEM29                1
    CNEH30    OBJ        2184.374512   DEM30                1
    CNEH31    OBJ        2184.374512   DEM31                1
    CNEH15    OBJ        2184.374268   DEM15                1
RHS
    RHS       DEM01              1.5   DEM02               69
    RHS       DEM03            520.5   DEM04       495.299805
    RHS       DEM05       754.099854   DEM06           1198.5
    RHS       DEM07      2188.799805   DEM08            594.4
    RHS       DEM09             51.4   DEM10             75.6
    RHS       DEM11       563.099854   DEM12      1533.099854
    RHS       DEM13       556.299805   DEM14            167.8
    RHS       DEM16       975.799805   DEM17           1642.9
    RHS       DEM18       347.099854   DEM19       262.799805
    RHS       DEM20               62   DEM21            395.9
    RHS       DEM22            691.7   DEM23            936.7
    RHS       DEM24       791.099854   DEM25            106.7
    RHS       DEM26            393.7   DEM27            173.7
    RHS       DEM28            466.2   DEM29               34
    RHS       DEM30            116.9   DEM31            103.4
    RHS       DEM15       552.099854
BOUNDS
 UP BOUNDS    D01DHT               1
 UP BOUNDS    D01GAZ               1
 UP BOUNDS    D02DHT               1
 UP BOUNDS    D02GAZ               1
 UP BOUNDS    D03DHT               1
 UP BOUNDS    D03GAZ               1
 UP BOUNDS    D04DHT               1
 UP BOUNDS    D04GAZ               1
 UP BOUNDS    D05DHT               1
 UP BOUNDS    D05GAZ               1
 UP BOUNDS    D06DHT               1
 UP BOUNDS    D06GAZ               1
 UP BOUNDS    D07DHT               1
 UP BOUNDS    D07GAZ               1
 UP BOUNDS    D08DHT               1
 UP BOUNDS    D08GAZ               1
 UP BOUNDS    D09DHT               1
 UP BOUNDS    D09GAZ               1
 UP BOUNDS    D10DHT               1
 UP BOUNDS    D10GAZ               1
 UP BOUNDS    D11DHT               1
 UP BOUNDS    D11GAZ               1
 UP BOUNDS    D12DHT               1
 UP BOUNDS    D12GAZ               1
 UP BOUNDS    D13DHT               1
 UP BOUNDS    D13GAZ               1
 UP BOUNDS    D14DHT               1
 UP BOUNDS    D14GAZ               1
 UP BOUNDS    D15DHT               1
 UP BOUNDS    D15GAZ               1
 UP BOUNDS    D16DHT               1
 UP BOUNDS    D16GAZ               1
 UP BOUNDS    D17DHT               1
 UP BOUNDS    D17GAZ               1
 UP BOUNDS    D18DHT               1
 UP BOUNDS    D18GAZ               1
 UP BOUNDS    D19DHT               1
 UP BOUNDS    D19GAZ               1
 UP BOUNDS    D20DHT               1
 UP BOUNDS    D20GAZ               1
 UP BOUNDS    D21DHT               1
 UP BOUNDS    D21GAZ               1
 UP BOUNDS    D22DHT               1
 UP BOUNDS    D22GAZ               1
 UP BOUNDS    D23DHT               1
 UP BOUNDS    D23GAZ               1
 UP BOUNDS    D24DHT               1
 UP BOUNDS    D24GAZ               1
 UP BOUNDS    D25DHT               1
 UP BOUNDS    D25GAZ               1
 UP BOUNDS    D26DHT               1
 UP BOUNDS    D26GAZ               1
 UP BOUNDS    D27DHT               1
 UP BOUNDS    D27GAZ               1
 UP BOUNDS    D28DHT               1
 UP BOUNDS    D28GAZ               1
 UP BOUNDS    D29DHT               1
 UP BOUNDS    D29GAZ               1
 UP BOUNDS    D30DHT               1
 UP BOUNDS    D30GAZ               1
 UP BOUNDS    D31DHT               1
 UP BOUNDS    D31GAZ               1
 UP BOUNDS    D32DHT               1
 UP BOUNDS    D32GAZ               1
 UP BOUNDS    D33DHT               1
 UP BOUNDS    D33GAZ               1
 UP BOUNDS    D34DHT               1
 UP BOUNDS    D34GAZ               1
 UP BOUNDS    D35DHT               1
 UP BOUNDS    D35GAZ               1
 UP BOUNDS    D36DHT               1
 UP BOUNDS    D36GAZ               1
 UP BOUNDS    D37DHT               1
 UP BOUNDS    D37GAZ               1
 UP BOUNDS    D38DHT               1
 UP BOUNDS    D38GAZ               1
 UP BOUNDS    D39DHT               1
 UP BOUNDS    D39GAZ               1
 UP BOUNDS    D40DHT               1
 UP BOUNDS    D40GAZ               1
 UP BOUNDS    D41DHT               1
 UP BOUNDS    D41GAZ               1
 UP BOUNDS    D42DHT               1
 UP BOUNDS    D42GAZ               1
 UP BOUNDS    D43DHT               1
 UP BOUNDS    D43GAZ               1
 UP BOUNDS    D44DHT               1
 UP BOUNDS    D44GAZ               1
 UP BOUNDS    D45DHT               1
 UP BOUNDS    D45GAZ               1
 UP BOUNDS    D46DHT               1
 UP BOUNDS    D46GAZ               1
 UP BOUNDS    D47DHT               1
 UP BOUNDS    D47GAZ               1
 UP BOUNDS    D48DHT               1
 UP BOUNDS    D48GAZ               1
 UP BOUNDS    D49DHT               1
 UP BOUNDS    D49GAZ               1
ENDATA
