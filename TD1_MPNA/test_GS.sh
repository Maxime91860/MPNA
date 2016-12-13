
rm -rf gram_schmidt.txt

for((NB_VECTEUR=2;$NB_VECTEUR<500;NB_VECTEUR=NB_VECTEUR+1))
	do
	echo "./gram_schmidt.pgr $NB_VECTEUR"
	./gram_schmidt.pgr $NB_VECTEUR
done