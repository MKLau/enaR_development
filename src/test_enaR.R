### enaR function output checking algorithm

rm(list=ls())
reset.base <- FALSE
if (!tail(strsplit(getwd(),split='/')[[1]],1) == 'enaR_development'){
    print('Working directory incorrect. Change to enaR_developemt')
}else{
    library(devtools)
    install_github('SEELab/enaR',ref='develop')
    library(enaR)
    data(oyster)
    test <- capture.output(
        Sys.time(),
        Sys.info(),
        sessionInfo(),
        as.extended(oyster),
        enaAll(oyster),
        enaMTI(oyster),
        environCentrality(as.matrix(oyster)),
        netOrder(oyster,6:1),
        read.wand('data/MDmar02_WAND.xls'),
        ssCheck(oyster),
        EcoNetWeb("Intertidal Oyster Reef Ecosystem Model"),
        enaAscendency(oyster),
        findPathLength(oyster),
        pack(unpack(oyster)$'F',unpack(oyster)$'z',unpack(oyster)$'r',unpack(oyster)$'e',unpack(oyster)$'y',unpack(oyster)$'X',unpack(oyster)$'living'),
        ShannonDiversity(unpack(oyster)$'X'),
        balance(oyster),
        enaControl(oyster),
        enaStorage(oyster),
        force.balance(oyster),
        read.EcoNet('data/pitcherCN.eco'),
        scc(oyster),
        unpack(oyster),
        TES(oyster),
        enaCycle(oyster),
        enaStructure(oyster),
        get.ns(oyster),
        read.enam('data/MDMAR02.xlsx'),
        write.EcoNet(oyster,'data/econet_test.txt'),
        TET(oyster),
        enaEnviron(oyster),
        enaTroAgg(oyster),
        get.orient(),
        read.nea('data/nea_test.txt'),
        set.orient('rc'),
        write.nea(oyster,'data/nea_test.txt'),
        as.bipartite(oyster,c(1,1,1,2,2,2)),
        eigenCentrality(as.matrix(oyster)),
        enaFlow(oyster),
        enaUtility(oyster),
        mExp(as.matrix(oyster)),
        read.scor('data/cone_spring.dat'),
        signs(as.matrix(oyster))
        )
    if (reset.base){
        fileout <- "data/test_base.txt"
    }else{
        fileout <- "data/test_new.txt"
    }
    fileConn <- file(fileout)
    writeLines(test, fileConn)
    close(fileConn)
    system(paste0('diff data/test_base.txt',' ','data/test_new.txt','> data/diff_base_new.txt'))
}
