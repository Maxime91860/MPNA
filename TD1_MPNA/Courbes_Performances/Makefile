

all :  scalabilite_forte.eps produit_scalaire PMV GS

produit_scalaire : produit_scal_croisement.eps produit_scal_temps_taille.eps

PMV : pmv_temps_taille.eps

GS : gram_schmidt.eps

%.eps : %.gplot
	gnuplot $<

clean :
	rm -rf *.eps
