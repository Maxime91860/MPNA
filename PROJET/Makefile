

all : code

exec :
	./CODE/bi_lanczos_seq.pgr 9 5

code : 
	make -C CODE	

clean_code : 
	make clean -C CODE
	rm log.txt

clean_output : 
	rm -rf OUTPUT/* 

clean_plot :
	rm -rf PLOT/* 

clean : clean_output clean_plot clean_code