#########################################################################################
### Helper Functions for RCPmod used in Foster, Hill & Lyons 2016 JRSS paper          ###
### Written by N. Hill July 2015.                                                     ###
#########################################################################################



#### poly_data----
# Save polynomial bases to transform prediction space

poly_data<-function(poly_vars,      #vector of predictor variable names to create orthogonal polynomials. Matched to column names in 'data'
                    degree,         #vector of same length as poly_vars specifying the polynomial degree for each predictor variable in poly_vars
                    id_vars,        #vector of ID variable names (Not transformed). Matched to column names in 'data'
                    sample_vars=NULL, #vector of sampling level variable names (e.g. Gear type). Matched to column names in 'data'
                    species_vars,    #vector of response variables names. Matched to column names in 'data'
                    offset=NULL,     #name of offset variable.  Matched to column names in 'data'
                    data, ...)
  {
  store_polys<-list()
  for(i in 1:length(poly_vars)){
    store_polys[[i]]<-poly(data[,poly_vars[i]], degree=degree[i])
    dimnames(store_polys[[i]])[[2]]<-paste0(poly_vars[i],seq(1:degree[i]))
  }
  names(store_polys)<-poly_vars

  rcp_data<-na.omit(cbind(subset(data, select= c(id_vars, sample_vars, offset, species_vars)),
                          do.call(cbind, store_polys)))

  return(list(rcp_data=rcp_data, poly_output=store_polys))
}


### poly_pred_space----
# Transforming prediction space using same basis as orthogonal polynomials used to build models
# Note: only accomodates one sampling term at the moment
# Note: offset isn't actually used in predict function that predicts RCP membership, but will keep as it might be useful to predict expected abundance of species at site.

poly_pred_space<-function(pred_space,             #dataframe containing variables
                          poly_output,            #extracted list of stored polynomial attribtutes from 'poly_data' function
                          offset_val=NULL,        #an offset value. Possibly mean of offset used in model building. Will be logged within function
                          offset_name=NULL,       #name of offset used in RCP model
                          sampling_vals=NULL,     #level of sampling factor for prediction
                          sampling_name=NULL,     #name of sampling factor used in RCP model
                          sampling_factor_levels=NULL )  #levels of sampling factor used in RCP model
  {

  # transform predictors using saved orthogonal polynomial attributes
  pred_polys<-list()
  vars<-names(poly_output)
  for( i in 1: length(vars)){
    pred_polys[[i]]<- predict( poly_output[[i]], pred_space[, names(pred_space) %in% vars[i]])
    dimnames(pred_polys[[i]])[[2]]<-dimnames(poly_output[[i]])[[2]]
  }
  pred_polys_df<-as.data.frame(do.call(cbind, pred_polys))

  #create offset term
  if(!is.null(offset_val)){
    pred_polys_df$offset<-log(offset_val)
    names(pred_polys_df)[ncol(pred_polys_df)]<-paste0("log(", offset_name, ")")
  }

  # create sampling variable.
  # only accommodates one sampling factor
  if(!is.null(sampling_vals)){
    reps<- length(sampling_vals)
    pred_polys_df$sampling<-factor(sampling_vals,levels=sampling_factor_levels)
    names(pred_polys_df)[ncol(pred_polys_df)]<-sampling_name
  }

  return(pred_polys_df)
}

