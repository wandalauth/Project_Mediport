#Variable Selection via Backwards Elimination 1000x

set.seed(1)
df_NoNa <- na.omit(df) #Only a full dataset can be used for backwards elimination

R <- 1000 #Number of runs
Variables <- matrix(nrow = R, ncol = 40) #Matrix to store final variables of backwards elimination
EPV <- rep(NA,R) #Matrix to store EPV per run


for (i in 1:R) {
  
  testdata_NAs <- df %>% 
    dplyr::filter(is.na(Age) | is.na(Residence) | is.na(`No of drugs at admission`) |
                    is.na(BMI) | is.na(ASA) | 
                    is.na(`hospital ward`) | is.na(`Hospitlisation last 12m`) |
                    is.na(CCI) | is.na(Complications) | is.na(Renal_function_group) |
                    is.na(Gender))
  
  #Splitting in train and test dataset
  traindata <- df_NoNa %>% dplyr::slice_sample(prop=0.70)
  testdata  <- dplyr::anti_join(df_NoNa, traindata, by = 'Patient No')
  
  EPV[i] <- sum(floor(traindata$Complications))
  
  testdata <- rbind(testdata, testdata_NAs)
  
  # Full model
  full.model <- glm(as.factor(Complications) ~., data = traindata[,-1], family = "binomial")
  
  # Stepwise regression model
  step.model <- stepAIC(full.model, direction = "backward", 
                        trace = F)
  
  #store final variables in matrix
  Variables[i,] <- c(names(step.model$coefficients)[-c(1)], rep(NA, ncol(Variables)-length(names(step.model$coefficients)[-c(1)])))
}


EPVMean <- data.frame(Mean = mean(EPV/12), SD = sd(EPV/12), Minimum = min(EPV/12), Maximum = max(EPV/12))

#Variables are only needed once, not for every level
Variables <- as.data.frame(Variables) %>%
  dplyr::mutate_all(funs(ifelse(.== "Renal_function_group30-59"  |
                                  . == "Renal_function_group60-89" |
                                  . == "Renal_function_group>89" , NA,.))) %>% 
  dplyr::mutate_all(funs(ifelse(.== "IntoleranceOther", NA,.))) %>% 
  dplyr::mutate_all(funs(ifelse(.== "IntoleranceNo", "Intolerance",.))) %>% 
  dplyr::mutate_all(funs(ifelse( . == "Renal_function_group15-29", "Renal_function_group",.))) %>% 
  dplyr::mutate_all(funs(ifelse( . == "`hospital ward`GY" |
                                   . == "`hospital ward`GF" |
                                   .== "`hospital ward`KG" |
                                   . == "`hospital ward`LA" | 
                                   . == "`hospital ward`OT" |
                                   . == "`hospital ward`SF" |
                                   . == "`hospital ward`UR" , NA,.))) %>% 
  dplyr::mutate_all(funs(ifelse(. == "`hospital ward`DE", "hospital ward",.))) %>%
  dplyr::mutate_all(funs(ifelse(. == "GenderW", "Gender",.))) %>% 
  dplyr::mutate(row = seq(1:R))


Variables <- as.matrix(Variables)

# Occurrence of total number of variables
N_Variables <- t(Variables) %>% as.data.frame() %>% 
  dplyr::summarise(across(everything(), ~sum(!is.na(.))-1))  %>% 
  t() %>% as.data.frame() %>%  group_by(V1) %>% 
  dplyr::summarise(N=n()) %>% 
  dplyr::mutate(RelFreq = round(N/R,3))

# Occurrence of variable combinations 
Variable_Combinations <-  Variables %>%  as.data.frame() %>%
  dplyr::mutate(row=NULL) %>% t() %>% as.data.frame() %>% 
  dplyr::summarise(across(everything(), ~ paste(.[!is.na(.)], collapse=',' ))) %>% 
  t() %>% as.data.frame() %>% 
  group_by(V1) %>% dplyr::summarise(N=n())


# Occurence of single variables
Occurrence_Variables <- as.data.frame(Variables) %>% 
  pivot_longer(-row, values_to = "Variable", names_to = "Names") %>% 
  filter(!is.na(Variable)) %>% 
  dplyr::mutate(Names = NULL) %>% 
  group_by(row, Variable) %>% 
  dplyr::summarise(N = n()) %>% 
  group_by(Variable) %>% 
  dplyr::summarise(N = n()) %>% 
  dplyr::mutate(RelFreq = round(N/R,5))


#Selection of the variables with an occurrence of >50% -> is also the most common combination
