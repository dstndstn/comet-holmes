all: holmes.pdf

quick:
	sed 's/\.eps/\.pdf/g' holmes.tex > tmp-holmes.tex
	pdflatex tmp-holmes.tex
	mv tmp-holmes.pdf holmes.pdf

arxiv: arXiv

arXiv: *.png *.jpg *.pdf *.tex
	tar cvf holmes_arXiv.tar holmes.tex \
		cbar.pdf \
		hparam-hist-0.pdf \
		hparam-hist-1.pdf \
		hparam-hist-2.pdf \
		xyzorbit-0.pdf \
		xyzorbit-1.pdf \
		xyzorbit-2.pdf \
		holmes-0004.jpg \
		holmes-0013.jpg \
		holmes-0044.jpg \
		holmes-0046.jpg \
		holmes-0052.jpg \
		holmes-0058-b.jpg \
		holmes-0059.jpg \
		holmes-0128.jpg \
		holmes-0270.jpg \
		holmes-0419.jpg \
		holmes-0559.jpg \
		holmes-0572-b.jpg \
		holmes-0618.jpg \
		holmes-0638.jpg \
		holmes-0639.jpg \
		holmes-0787.jpg \
		holmes-0833.jpg \
		holmes-1219-b.jpg \
		holmes-1300.jpg \
		holmes-1632-b.jpg \
		holmes-1637.jpg \
		holmes-1641.jpg \
		holmes-2050.jpg \
		holmes-coadd-1.jpg \
		holmes-coadd-2.pdf \
		holmes-footprint-1.png \
		holmes-footprint-2.pdf \
		exif1.pdf \
		exif2.pdf \
		exif4.pdf \
		exif6.pdf \
		exif8.pdf \
		exif10.pdf \
		traj-default-nw064-bw008-emp-3000.pdf \
		hya-coadd-1.jpg \
		hya-coadd-2.pdf \
		hya-footprint-1.png \
		hya-footprint-2.pdf \
		param-hist-0.pdf \
		param-hist-1.pdf \
		param-hist-2.pdf \
		param-hist-3.pdf \
		param-hist-4.pdf \
		param-hist-5.pdf \
		lnp.pdf \
		pexif.pdf \
		pgood.pdf \
		eta.pdf \
		holmes-wco-1.png \
		hya-wco-1.png
	gzip --best -fv holmes_arXiv.tar

%.pdf: %.tex
	sed 's/\.eps/\.pdf/g' $< > tmp-$<
	pdflatex tmp-$<
	bash -c " ( grep Rerun tmp-$*.log && pdflatex tmp-$< ) || echo noRerun "
	bash -c " ( grep Rerun tmp-$*.log && pdflatex tmp-$< ) || echo noRerun "
	mv tmp-$@ $@
