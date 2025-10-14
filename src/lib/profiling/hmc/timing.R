args = commandArgs(trailingOnly=TRUE)

mnl <- function () {
  mnl <- list()
  class(mnl) <- append(class(mnl), 'mnl')
  return (mnl)
}

create_mnl<-function(.mnl=mnl(),name,type){
  stopifnot(inherits(.mnl, 'mnl'))
  
  .mnl$name<-name
  .mnl$type<-type
  
  .mnl$sub_func<- list()
  
  # .mnl$hb_name<-paste0()
  # .mnl$deriv_name<-c()
  # .mnl$hb_name<-c()
  # 
  
  .mnl$hb_time<-c()
  .mnl$deriv_time<-c()
  .mnl$hb_time<-c()
  # class(.mnl) <- append(class(.mnl), 'mnl_meta')
  return(.mnl)
} 

sub_func <- function () {
  sub_func <- list()
  class(sub_func) <- append(class(sub_func), 'sub_func')
  return (sub_func)
}

sub_func_create<-function(.sub_func=sub_func(),name,time_in , level ){
  stopifnot(inherits(.sub_func, 'sub_func'))
  
  .sub_func$name<-name
  .sub_func$time<-time_in
  .sub_func$level<-level
  # class(.mnl) <- append(class(.mnl), 'mnl_meta')
  return(.sub_func)
} 



string<-paste("grep -E \"initialising monomial with type|monomial named\" ",args[1]," > mnl_tmp.txt" )
system(string) 
file <- "mnl_tmp.txt"
infile<- read.table(file, fill = TRUE ,comment.char = "")

df<-data.frame("name"=c(),"type"=c())
mnl_type<- c()
mnl_name<- c()
for (i  in c(1:length(infile[,1]))  ){
  type<-infile[i,5]
  name<-infile[i,3]
  #cat(i,name, type,"\n\n")
  if (infile[i,4]=="type" ){
    mnl_type<-c(mnl_type ,type )
    if (type=="GAUGE"){
      mnl_name<-c( mnl_name, "GAUGE" )
    }
  }
  else if (infile[i,2]=="named" )
    mnl_name<-c(mnl_name ,name )
  
}
df<-data.frame("name"=mnl_name,"type"=mnl_type)
mnl_list<- vector("list",length(mnl_name))
mnl_list<-list()
for (i  in c(1:length(mnl_name))  ){
  #mnl_list<- append(mnl_list,create_mnl( name=mnl_name[i], type=mnl_type[i]))
  mnl_list[[i]]<- create_mnl( name=mnl_name[i], type=mnl_type[i])
}

acc_names<-paste0()
deriv_names<-c()
hb_names<-c()
# df$hb_times<-c(0)
# df$deriv_times<-c(0)
# df$acc_times<-c(0)




#grep the lines with time output 
# tac revese the lines in the file
string<-paste0("grep  \"Time\" ",args[1]," | tac > tmp" )########################################uncomment here and line below
system(string) 
file <- "tmp"
infile<- read.table(file, fill = TRUE ,comment.char = "")

time_hb<-c(list())
time_deriv<-c(list())
time_acc<-c(list())

index_mon <- function(name){
  for(i in c(1:length(mnl_list))){
    if (mnl_list[[i]]$name == name)
      return(i)
  }
  return(-1)  
}

current_level<-1e+6
current_id_mon<-1
current_mnl<-"non_existing_mnl"

current_func<-c("non_existing_func")
current_id_fun<-c(1)
nesting<-1


add_sub_func <- function(mnl_sub_func ,func_tree, secs, level){
  
  
  Nfcheck <- length(mnl_sub_func)
  Nfin <- length(func_tree)
  #if no function is found
  if ( Nfcheck <= 0 ){
    #tmp <- list()
    #cat("empty sub_func list create:",func_tree[1],"\n")
    tmp <- list(sub_func_create(name=func_tree[1],time_in = secs,level = level))
    return(tmp)
  }
  #if we are at the deeepest level
  else if ( Nfin == 1 )  {
    found=FALSE
    for (j in c(1:Nfcheck)){
      if(mnl_sub_func[[j]]$name == func_tree[1] ){
        tmp <- mnl_sub_func
        tmp[[j]]$time <- c(tmp[[j]]$time, secs)
        #cat("found existing function:",func_tree[1],"\n")
        return(tmp)
        found=TRUE
      }
    }
    if(!found){
      tmp <- mnl_sub_func
      #print("creating new")
      tmp[[Nfcheck+1]] <- sub_func_create(name=func_tree[1],time_in=secs,level=level)
      return(tmp)
    }
  }
  # resursion in level is not the deepest
  else if(Nfin>1){
    found=FALSE
    for (j in c(1:Nfcheck)){
      if(mnl_sub_func[[j]]$name == func_tree[1] ){
        tmp <- mnl_sub_func
        tmp[[j]]$sub_func <- add_sub_func(tmp[[j]]$sub_func ,func_tree[-1], secs, level)
        return(tmp)
        found=TRUE
      }
    }
    stopifnot(!found);
    
    
  }
}


N <- length(infile[,1])
init=TRUE
for(i in c((1):N)){######################################################################### debug 1:N
  
  name <- infile[i,2]
  name <- gsub(":","",name)
  id <- index_mon(name)
  
  if (id == -1) {
    
    current_level<-1e+6
    current_id_mon<-1
    current_mnl<-"non_existing_mnl"
    
    current_func<-c("non_existing_func")
    current_id_fun<-c(1)
    next;
  }
  func_name <- infile[i,6]
  sec <- as.numeric(infile[i,7])
  level <- as.integer(infile[i,10])
  #cat(name ,"   ",func_name, "   ",sec, "   ",level ,"\n" )
  if (init ){
    current_level <- level
    current_id_mon <- id
    current_mnl <- name
    
    current_func <- c(func_name)
    current_id_fun <- c(1)
    init=FALSE
  }
  else{
      if(current_mnl != name){
        #reset current_func
        #cat("changing monomial to id=",id, "\n")
        current_func <- c(func_name)
      }
      else{
      
        if(level>current_level){
          #cat("moving from sublevel",current_level,"to",level, "\n")
          
          current_func<-c(current_func,func_name)
          
        }
        else if(level<current_level){
          #cat("moving from sublevel",current_level,"to",level, "\n")
          if (length(current_func) < 2){
            current_func <- c(func_name)
          }
          else{ 
            # remove the last up to the new level
            
            ln <- c(1:(length(current_func) - current_level + level-1 ) )
            current_func<-c(current_func[ln] ,func_name )
            #current_id_fun<-c(current_id_fun(ln),func_name)
          }
        }
        else if(level==current_level){
          #cat("adding func to the same level sublevel",current_level, "\n")
          
          ln<-length(current_func)
          current_func[ln]<- func_name
        }  
      }
  
        
        
  }
  current_level<-level
  current_id_mon<-id
  current_mnl<-name
  mnl_list[[id]]$sub_func <- add_sub_func(mnl_list[[id]]$sub_func, current_func, sec, level)
  
  
}
save(mnl_list, file = "profile.RData")

rmarkdown::render("profile.Rmd")