require(bio3d)
require(rgl)
# require(Rpdb)
require(ggplot2)


samples_name <- list.files(path = "experimental-pdb/")


pdbs_exp <- pdbaln(paste("experimental-pdb/", samples_name, sep = ""))
pdbs_exp


a1_st1 <- read.pdb("experimental-pdb/URA3_a1_st1.B99990004.pdb")


gdk <- bio3d::read.pdb("tmp/3gdk.pdb")


pdbsplit(get.pdb("3GDK", URLonly=T))


gdk_a <- read.pdb("split_chain/3GDK_A.pdb")

gdk_a$sheet

pdbseq(gdk_a) %>% length()


