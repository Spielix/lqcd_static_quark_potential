reset

set encoding iso_8859_1
#set grid
set samples 10000
#set pointsize 0.1
set format '$%g$'
#set xtics ('0' 0, '' 1000, '2000' 2000, '' 3000, '4000' 4000, '' 5000, '6000' 6000, '' 7000, '8000' 8000)
set terminal epslatex size 15.0cm, 10.0cm color colortext
set fit results
#set fit errorvariables
set pointsize

set autoscale
set xlabel '$\frac{\beta - \beta_c}{\beta_c}N_\sigma^{1/\nu}$'
set ylabel '$\frac{\langle L^4 \rangle}{\langle L^2 \rangle^2}$'
set label 2 at -0.25, 1.3 '$\beta_c=\num{2.2977+-0.0002}$, $\nu=\num{0.64+-0.02}$, reduced $\chi^2 = \num{1.05}$' front

A = 0.5981421; B = -1.1461195; C = 1.6571032; beta_c = 2.2977317; nu = 0.6354707; dbeta_c = 0.0001333353;
dnu = 0.0131096580; cor_bcnu = -0.2355493;
h(x) = A*x**2 + B*x + C
#fit [40:130] h(x) './kalibration_4.csv' using ($1-7.12):2:(sqrt($2)) zerror via I,J,K

set output './fin_size_scaling.tex'
plot 'fin_size_scaling_L08.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 -beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 -($1 - beta_c) * $1 * $3**(2/nu-1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2 * $16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma =\  8$' with xyerrorbars,\
    'fin_size_scaling_L09.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 -beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2 *$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma =\  9$' with xyerrorbars,\
    'fin_size_scaling_L10.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 10$' with xyerrorbars,\
    'fin_size_scaling_L11.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 11$' with xyerrorbars,\
    'fin_size_scaling_L12.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 12$' with xyerrorbars,\
    'fin_size_scaling_L13.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 13$' with xyerrorbars,\
    'fin_size_scaling_L14.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 14$' with xyerrorbars,\
    'fin_size_scaling_L15.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 15$' with xyerrorbars,\
    'fin_size_scaling_L16.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 16$' with xyerrorbars,\
    h(x) title 'Polynom zweiter Ordnung'

set output "beta_scan_order.tex"
set xlabel "$\beta$"
set ylabel "$\langle |L| \rangle$"
plot "fin_size_scaling_L16.csv" using 1:8:9 notitle with yerrorbars

set output "beta_scan_order_ctime.tex"
set ylabel "$\tau_{\text{int},|L|}$"
plot "fin_size_scaling_L16.csv" using 1:10 notitle with linespoints
