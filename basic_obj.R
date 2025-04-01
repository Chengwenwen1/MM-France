library(ggplot2)
mytheme <- theme_bw() +
  theme(#legend.position = "right",
    axis.text.x=element_text(angle=0,hjust=0.5),
    text = element_text(size=12,face="bold"),
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(), 
    strip.text = element_text(size=12,face="bold"),
    strip.background = element_blank(),
    axis.line = element_line(color = 'black'))
library(RColorBrewer)
type.col <- colorRampPalette(c(brewer.pal(n=3,name = 'Set3'),
                               brewer.pal(n=8,name='Dark2')))(19)
colorRampPalette()(15)


big.marker <- c('Cd34','Avp','Prss57',#HSC
                'Klf1','Gypa','Slc4a1',#'CD36','TSPO2',#EPC
                'Hba1', 'Hbb', #Erythrocytes
                'Cxcl12','Col6a1','Col1a1','Col3a1',#Stromal cell
                'Cd3d','Cd3e','Cd4','Cd40lg','Sell',#naive
                'Ccr7','Lef1',#naive
                'Cd8a','Cd8b',#CD8+ T cells
                'Gzmk',#'GZMH',#CD8 mem
                'Gnly','Gzmb',#CD8 E
                'Klrf1','Klrc1','Ncam1',#NK cell#gnly'KLRB1',
                'Znf683',#NKT
                'Vpreb1','Tcl1a','Cd19','Cd79a','Ms4a1','Tcl1a',#B cells
                'Lyz', #Myeloid cells
                'Cd14', #classic Monocytes, 
                #'SEMA6B',
                'Fcgr3a', #non classic Monocytes
                'Retn','Mki67',#prMono
                #'ELANE',#gmp
                'Clec9a',##cDC1 
                "Clec10a",'Fcer1a','Cd1c',##cDC2
                "Il3ra","Lilra4",#pDC
                'Sdc1','Ighg1','Mzb1'#Plasma cells
)


