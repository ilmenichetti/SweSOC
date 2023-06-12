# load ggvenn package
library("ggvenn")

# use data frame as input
M <-tibble(value=c(1,3,2,7,5),'TEA'=c(TRUE,FALSE,TRUE,FALSE,FALSE),
           'COFFEE'=c(TRUE,TRUE,FALSE,FALSE,TRUE),
           'JUICE'=c(TRUE,FALSE,FALSE,TRUE,TRUE))

# create Venn diagram and display all sets
ggvenn(M)


model_names<-c("Olson","ICBM", "Century", "Yasso", "Yasso07", "RothC", "Porous", "3PFeedback")

venn_models<-mat.or.vec(length(model_names), 4)
rownames(venn_models)<-model_names
colnames(venn_models)<-c("Serial", "Parallel", "Inert", "Feedback")

venn_models[1,]=c(1,0,0,0)
venn_models[2,]=c(1,0,0,0)
venn_models[3,]=c(0,1,0,0)
venn_models[4,]=c(1,0,0,0)
venn_models[5,]=c(1,0,0,0)
venn_models[6,]=c(1,0,0,1)

# use data frame as input

venn_models_tib <-tibble(value=c("Serial", "Parallel", "Inert", "Feedback"),
           "Olson"=c(T,F,F,F),
           "ICBM"=c(T,F,F,F),
           "Century"=c(T,F,F,F),
           "Yasso"=c(T,F,F,F),
           "Yasso07"=c(T,F,F,F),
           "RothC"=c(T,F,F,F),
           "Porous"=c(T,F,F,F),
           "3PFeedback"=c(T,F,F,F))


venn_models_tib <-tibble(value=model_names,
                         "Serial"=c(T,T,F,T,T,T,T,T),
                         "Parallel"=c(F,T,T,F,F,F,F,F),
                         "Inert"=c(F,F,F,F,F,F,T,T),
                         "Feedback"=c(F,F,T,F,F,F,F,T))

# create Venn diagram and display all sets
ggvenn(venn_models_tib)



x<- list()
x$A <- as.character(c("Per_36","Cent","CeM","vDG","LAVL","RSGd"))
x$B <- as.character(c("vCA1","DLE","Per_36","vDG","DIE","Per_35"))
x$C <- as.character(c("vCA1","Cg1","LAVL", "RSGc", "RSGd","Per_35","Per_36"))
x$D <- as.character(c("Por","Cg1","RSGc","LAVL","Per_35","RSGd","Per_36"))


require(VennDiagram)
v0 <-venn.diagram(x, lwd = 3, col = c("red", "green", "orange", "blue"),
                  fill = c("red", "blue", "green", "orange"), apha = 0.5, filename = NULL)
grid.draw(v0)
overlaps <- calculate.overlap(x)
overlaps <- rev(overlaps)
for (i in 1:length(overlaps)){
  v0[[i+8]]$label <- paste(overlaps[[i]], collapse = "\n")
}

grid.newpage()
grid.draw(v0)
