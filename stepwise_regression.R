##### test stepwise regression model

regtest <- read.csv("C:/Users/Meta/Desktop/Data for R/regtest.csv")

regtest1 <- step(lm(AMGO~mowtime+ACRU+AGPA+APCA,data=regtest),direction="both")

####check for redundancy
cor(vegspecies_reg, method="spearman")

summary(regtest1)

######### stepwise regression model for 2015 by bird species##########

datastepreg2015 <- read.csv("C:/Users/Meta/Desktop/Data for R/datastepreg2015.csv")

###AMGO
AMGOtest <- step(lm(AMGO ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(AMGOtest)

###AMRO
AMROtest <- step(lm(AMRO ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(AMROtest)

###BAOR
BAORtest <- step(lm(BAOR ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(BAORtest)

###BARS
BARStest <- step(lm(BARS ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(BARStest)

###BAWW
BAWWtest <- step(lm(BAWW ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(BAWWtest)

###BCCH
BCCHtest <- step(lm(BCCH ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(BCCHtest)

###BGGN
BGGNtest <- step(lm(BGGN ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(BGGNtest)

###BHCO
BHCOtest <- step(lm(BHCO ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(BHCOtest)

###BLJA
BLJAtest <- step(lm(BLJA ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(BLJAtest)

###BRTH
BRTHtest <- step(lm(BRTH ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(BRTHtest)

###BWWA
BWWAtest <- step(lm(BWWA ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(BWWAtest)

###CEDW
CEDWtest <- step(lm(CEDW ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(CEDWtest)

###COGR
COGRtest <- step(lm(COGR ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(COGRtest)

###COYE
COYEtest <- step(lm(COYE ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt, data=datastepreg2015),direction="both")

summary(COYEtest)

###EABL
EABLtest <- step(lm(EABL ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(EABLtest)

###EAKI
EAKItest <- step(lm(EAKI ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(EAKItest)

###EATO
EATOtest <- step(lm(EATO ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(EATOtest)

###FISP
FISPtest <- step(lm(FISP ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(FISPtest)

###GCFL
GCFLtest <- step(lm(GCFL ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(GCFLtest)

###GRCA
GRCAtest <- step(lm(GRCA ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(GRCAtest)

###HOWR
HOWRtest <- step(lm(HOWR ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(HOWRtest)

###INBU
INBUtest <- step(lm(INBU ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(INBUtest)

###MODO
MODOtest <- step(lm(MODO ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(MODOtest)

###NOCA
NOCAtest <- step(lm(NOCA ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(NOCAtest)

###NOFL
NOFLtest <- step(lm(NOFL ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(NOFLtest)

###OROR
ORORtest <- step(lm(OROR ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(ORORtest)

###OVEN
OVENtest <- step(lm(OVEN ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(OVENtest)

###PUMA
PUMAtest <- step(lm(PUMA ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(PUMAtest)

###RBGR
RBGRtest <- step(lm(RBGR ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(RBGRtest)

###RBWO
RBWOtest <- step(lm(RBWO ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(RBWOtest)

###RWBL
RWBLtest <- step(lm(RWBL ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(RWBLtest)

###SOSP
SOSPtest <- step(lm(SOSP ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(SOSPtest)

###SWSP
SWSPtest <- step(lm(SWSP ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(SWSPtest)

###TRES
TREStest <- step(lm(TRES ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(TREStest)

###TUTI
TUTItest <- step(lm(TUTI ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(TUTItest)

###VEER
VEERtest <- step(lm(VEER ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(VEERtest)

###WAVI
WAVItest <- step(lm(WAVI ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(WAVItest)

###WBNU
WBNUtest <- step(lm(WBNU ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(WBNUtest)

###WEVI
WEVItest <- step(lm(WEVI ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(WEVItest)

###WIFL
WIFLtest <- step(lm(WIFL ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(WIFLtest)

###WITU
WITUtest <- step(lm(WITU ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(WITUtest)

###WOTH
WOTHtest <- step(lm(WOTH ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(WOTHtest)

###YBCU
YBCUtest <- step(lm(YBCU ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(YBCUtest)

###YEWA
YEWAtest <- step(lm(YEWA ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(YEWAtest)

###YTVI
YTVItest <- step(lm(YTVI ~ ACRU + APCA + ARVU + CARE + COAM + CORA + DICL + EUGR + GAMO + JUNC + LIST + LOMO + LYLI + LYSA + MACA+ MIVI + MYPE + ONSE + PAQU + PHAR + POSI + PYVI + QUPA + ROMU + ROPA + RUBU + SCSC+ SCCY + SOLI + SPLA + SPTO + TORA + VACC + VENO + VIDE + VILA + mowtime + GCVR.Grass + GCVR.Forb + GCVR.ShrubTree + GCVR.Litter + GCVR.Wood + GCVR.Moss + GCVR.Water + VegHgt.Avg + VegHgt.Stdv + LD.Avg + LD.Stdv + AvgCanopy + ShrubSTM + TreeSTM + Avg.TreeDist + Avg.TreeHgt + Avg.ShrubDist + Avg.ShrubHgt + Avg.ShrubWdt,data=datastepreg2015),direction="both")

summary(YTVItest)



