cat("\014") #Clears the console
rm(list=ls()) #Clears the environment

library("rmarkdown") #loads rmarkdown
#Set your working directory to wherever your data is stored
setwd("Your Working Directory")


#BE SURE YOU HAVE MADE AN EMPTY FOLDER NAMED "Assigned Formulas" IN THE FOLDER WITH YOUR DATA
#IF YOU HAVE MULTIPLE SUB-FOLDERS CONTAINING .CSV FILES WITHIN YOUR DATA FOLDER, THIS FUNCTION WILL PULL ALL OF THEM
#SO BE CAREFUL IF YOU DON'T WANT TO TRY PROCESSING ALL THE FOLDERS.

#YOU WILL NEED A COMMON EXTENSION AT THE END OF ALL THE FILES, BY DEFAULT IT IS HAVE _MS IT CAN BE CHANGED IF DESIRED,
#BUT MAY CAUSE PROBLEMS IF NOT DONE CORRECTLY.
folder=paste0(getwd(),"/")    #Pulls your working directory
file_list=list.files(path=folder, pattern="*_MS.csv")


#MAKE SURE THE MARKDOWN IN THE "RENDER" CALL IS THE ONE YOU ARE PLANNING ON USING
#ALSO ENSURE THE MARKDOWN IS IN THE FOLDER WITH YOUR DATA

for (i in 1:length(file_list)){
  render("Formula Assignment Markdown.rmd",output_file = paste0(substr(file_list[i],1,nchar(file_list[i])-4),".html"),
         params=list(data=file_list[i]))
}
