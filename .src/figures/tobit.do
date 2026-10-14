cd "d:\Yandex.Disk\03. Преподавание\Эконометрика\.fig"
infile x y using "tobit.txt", clear

gen yc = 0
replace yc = y if y > 0

matrix COEFS = J(1, 2, .)

reg y x
mat COEFS = COEFS \ e(b)

reg yc x
mat COEFS = COEFS \ e(b)

reg y x if y > 0
mat COEFS = COEFS \ e(b)

tobit y x, ll(0)
mat COEFS = COEFS \ e(b)[1, 1..2]

mat COEFS = COEFS[2..5, 1..2]

clear
svmat COEFS, names(col)
export delimited using "tobit_betas.txt", replace novar delimiter(tab)
