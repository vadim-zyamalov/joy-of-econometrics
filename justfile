default := '/d/Yandex.Disk/03.\ Преподавание/Эконометрика'

sync path=default:
    ./sync.sh {{path}} .
    ./sync.sh . {{path}}

sync-bib:
    cp .data/bibliography.bib МГТУ\ им.\ Баумана/.src/bibliography.bib
    cp .data/bibliography.bib МГУ/.src/bibliography.bib
    cp .data/bibliography.bib МФТИ/.src/bibliography.bib
    cp .data/bibliography.bib .lehrbuch/Межобъектные\ данные/bibliography.bib
    cp .data/bibliography.bib .lehrbuch/Временные\ ряды/bibliography.bib