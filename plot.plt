reset

set encoding iso_8859_1
#set grid
set samples 10000
#set pointsize 0.1
set format '$%g$'
#set xtics ('0' 0, '' 1000, '2000' 2000, '' 3000, '4000' 4000, '' 5000, '6000' 6000, '' 7000, '8000' 8000)
set terminal epslatex size 15.0cm, 9.0cm color colortext
set fit results
#set fit errorvariables
set pointsize

set autoscale
set xlabel '$\frac{\beta - \beta_c}{\beta_c}N_\sigma^{1/\nu}$'
set ylabel '$\frac{\langle L^4 \rangle}{\langle L^2 \rangle^2}$'
set label 2 at -0.25, 1.3 '$\beta_c=\num{2.2977+-0.0002}\qquad\nu=\num{0.64+-0.02}$' front
set label 3 at -0.25, 1.4 'reduziertes $\chi^2 = \num{1.05}$' front

A = 0.5918199; B = -1.1450902; C = 1.6570540; beta_c = 2.2977362; nu = 0.635275; dbeta_c = 0.0001339531;
dnu = 0.0132487792; cor_bcnu = -0.2198784;
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
    h(x) title 'angepasste Parabel'


set label 2 at -0.15, 1.3 '$\beta_c=\num{2.2983+-0.0002}\qquad\nu=\num{0.72+-0.02}$' front
set label 3 at -0.15, 1.4 'reduziertes $\chi^2 = \num{2.53}$' front
B = -1.6625897; C = 1.6515988; beta_c = 2.2983120; nu = 0.7199421; dbeta_c = 0.0001510794;
dnu = 0.0151738785; cor_bcnu = -0.5583943;
h1(x) = B*x + C
#fit [40:130] h(x) './kalibration_4.csv' using ($1-7.12):2:(sqrt($2)) zerror via I,J,K

set output './fin_size_scaling_straight.tex'
plot 'fin_size_scaling_L08.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 -beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 -($1 - beta_c) * $1 * $3**(2/nu-1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2 * $16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma =\  8$' with xyerrorbars,\
    'fin_size_scaling_L09.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 -beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2 *$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma =\  9$' with xyerrorbars,\
    'fin_size_scaling_L10.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 10$' with xyerrorbars,\
    'fin_size_scaling_L11.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 11$' with xyerrorbars,\
    'fin_size_scaling_L12.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 12$' with xyerrorbars,\
    'fin_size_scaling_L13.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 13$' with xyerrorbars,\
    'fin_size_scaling_L14.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 14$' with xyerrorbars,\
    'fin_size_scaling_L15.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 15$' with xyerrorbars,\
    'fin_size_scaling_L16.csv' using ((($1 - beta_c) / beta_c) *$3**(1 / nu)):($16 / $12**2):(sqrt((($1 - beta_c) * dnu * $3**(1 / nu - 1) / nu / beta_c)**2 + ($1 * dbeta_c * $3**(1 / nu) / beta_c**2)**2 - ($1 - beta_c) * $1 * $3**(2 / nu - 1) / beta_c**3 / nu * cor_bcnu * dnu * dbeta_c)):(sqrt(($17 / $12**2)**2 + (2*$16 * $13 / $12**3)**2 - 4 * $20 * $16 * $17 * $13 / ($12**5 * $15 * $19))) title '$N_\sigma = 16$' with xyerrorbars,\
    h1(x) title 'angepasste Gerade'




set output "overr_metropolis_comp.tex"
set xlabel '$\beta$'
set ylabel '$\langle |L| \rangle$'
set xrange [2.289:2.311]
set yrange [0.09:0.17]
unset label 2
C1=-5.611339;B1=2.491698;C2=-5.664823;B2=2.514981;
f1(x)=C1+B1*x
f2(x)=C2+B2*x
plot "fin_size_scaling_L12_pure_metro_new.csv" using 1:8:9 title 'Nur Metropolis Updates: Daten' with yerrorbars,\
	"fin_size_scaling_L12.csv" using 1:8:9 title 'Metropolis und Overrelaxation Updates: Daten' with yerrorbars,\
	f1(x) title 'Nur Metropolis Updates: Anpassung $\chi^2=\num{1.26}$',\
	f2(x) title 'Metropolis und Overrelaxation Updates: Anpassung $\chi^2=\num{1.85}$'
set autoscale

set output "temperature_scan_order.tex"
set xlabel '$N_\tau$'
set xrange [1.5:8.5]
plot "temperature_scan.csv" using 2:8:9 notitle with yerrorbars

set output "temperature_scan_order_ctime.tex"
set ylabel '$\tau_{\text{int},|L|}$'
set yrange [0:40]
plot "temperature_scan.csv" using 2:10:21 notitle with yerrorbars
set autoscale

set output "beta_scan_order.tex"
set xlabel '$\beta$'
set ylabel '$\langle |L| \rangle$'
set xrange [1.95:2.55]
plot "beta_scan.csv" using 1:8:9 notitle with yerrorbars

set output "beta_scan_order_ctime.tex"
set ylabel '$\tau_{\text{int},|L|}$'
set yrange [0:45]
plot "beta_scan.csv" using 1:10:21 notitle with yerrorbars
