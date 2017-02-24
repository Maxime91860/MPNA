

all : code

exec :
	./CODE/bi_lanczos_seq.pgr 9 5

test :
	./SCRIPTS_TEST/tests_taille.sh

code : 
	make -C CODE	

plot :
	gnuplot SCRIPTS_VISU/*

visu :
	evince PLOT/*

clean_code : 
	make clean -C CODE
	rm log.txt

clean_output : 
	rm -rf OUTPUT/* 

clean_plot :
	rm -rf PLOT/* 

clean : clean_code