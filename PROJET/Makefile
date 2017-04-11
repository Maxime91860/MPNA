

all : code

code : 
	make -C CODE	

exemple : code
	./CODE/bi_lanczos_seq.pgr MATRICES/bp_400.mtx 10

clean : 
	make clean -C CODE



# Pour le tracer des courbes
test :
	./SCRIPTS_TEST/tests_taille.sh

plot :
	gnuplot SCRIPTS_VISU/*

visu :
	evince PLOT/*

clean_plot :
	rm -rf PLOT/* 

clean_output : 
	rm -rf OUTPUT/* 