all:

	make sim
	make convert_to_formats

sim:

	oommf boxsi test_mif.mif

OUTF = "omfs_formats"
OMFFILE = "test-Oxs_MinDriver-Magnetization-00-0000000.omf"
convert_to_formats:
	mkdir -p $(OUTF)
	oommf avf2ovf -format text -grid irreg $(OMFFILE) $(OUTF)/test_irreg_txt.omf; \
	oommf avf2ovf -format text -grid reg   $(OMFFILE) $(OUTF)/test_rect_txt.omf; \
	oommf avf2ovf -format b4 -grid irreg   $(OMFFILE) $(OUTF)/test_irreg_b4.omf; \
	oommf avf2ovf -format b4 -grid reg     $(OMFFILE) $(OUTF)/test_rect_b4.omf; \
	oommf avf2ovf -format b8 -grid irreg   $(OMFFILE) $(OUTF)/test_irreg_b8.omf; \
	oommf avf2ovf -format b8 -grid reg     $(OMFFILE) $(OUTF)/test_rect_b8.omf

clean:

	rm -f test_*.omf

	rm -f test.odt
