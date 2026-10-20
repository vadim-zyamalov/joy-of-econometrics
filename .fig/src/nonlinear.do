cd "d:\Yandex.Disk\03. Преподавание\Эконометрика\.fig"
infile x y using "nonlinear.txt", clear

gen xsq = x^2
reg y x xsq

clear
svmat e(b), names("x")
export delimited using "nonlinear_betas.txt", replace novar delimiter(tab)
