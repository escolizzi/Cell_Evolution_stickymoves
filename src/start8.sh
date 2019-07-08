rm -r data_cellcount*.txt data_film2*/ backup*/;

./cell_evolution ../data/cell_evolution.par &
./cell_evolution ../data/cell_evolution.par -datafile data_cellcount2.txt -backupdir backup2/ -datadir data_film22 &
./cell_evolution ../data/cell_evolution.par -datafile data_cellcount3.txt -backupdir backup3/ -datadir data_film23 &
./cell_evolution ../data/cell_evolution.par -datafile data_cellcount4.txt -backupdir backup4/ -datadir data_film24 &
./cell_evolution ../data/cell_evolution.par -datafile data_cellcount5.txt -backupdir backup5/ -datadir data_film25 &
./cell_evolution ../data/cell_evolution.par -datafile data_cellcount6.txt -backupdir backup6/ -datadir data_film26 &
./cell_evolution ../data/cell_evolution.par -datafile data_cellcount7.txt -backupdir backup7/ -datadir data_film27 &
./cell_evolution ../data/cell_evolution.par -datafile data_cellcount8.txt -backupdir backup8/ -datadir data_film28 &

