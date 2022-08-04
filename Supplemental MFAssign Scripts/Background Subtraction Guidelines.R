###Background Subtraction After Formula Assignment

#For the purpose of this example, "Blank" corresponds to the assigned mass list for the blank sample 
#and "Sample" corresponds to the assigned mass list for the sample (Both are outputs from MFAssign)

#First I would recommend selecting only a few columns from the blank and 
#renaming a couple columns to make alignment easier

Blank2 <- Blank2 %>% select(formula, exp_mass, abundance) %>% 
  rename(blank_mass = exp_mass, blank_abund = abundance)

#Next you can align this Blank subset to the actual sample dataframe

Combo <- left_join(Sample, Blank2, by = "formula")

#Now you can make a decision about how you want to remove peaks from the sample.

#You can do it with presence, so if a peak is in the blank at all you remove it from the sample

Final <- Combo %>% filter(is.na(blank_mass)) #This removes all formulas that matched between blank and sample

#Sometimes a peak is very small in the blank but larger in the sample so it may not be reasonable to 
#remove the peak just because it matches. In this case you can do removal by ratio (for example it
#needs to be 10x larger in the sample than in the blank to be kept)

Final <- Combo %>% mutate(Ratio = abundance/blank_abund) %>% filter(Ratio > 10)

#You can also just subtract the blank abundance and remove those where the new abundance goes below 
#zero. The way to do it is the same as for the ratio method, just with a different equation.

###############################################
################################################
#The other way to do it is before assignment. To do this I have typically ended up using the 
#"expand.grid" function in R which matches all values between two dataframes and from there you can 
#filter it down to what you want to see.

library(MFAssignR)
library(dplyr)
Sample <- Raw_Neg_ML %>% filter(intensity > 1800) %>% rename(Samp_mass = m.z)
Blank <- Raw_Neg_ML %>% filter(intensity > 1800 & intensity < 2100) %>% rename(Blank_mass = m.z)
#This creates the big expanded data frame that matches all masses
Combo <- expand.grid(Sample$Samp_mass, Blank$Blank_mass)

#This renames the colums so it is easier to keep track of them.
Combo <- Combo %>% rename(Samp_mass = Var1, Blank_mass = Var2)

#Now we need to add the abundances for the sample and blank back in to help with deciding what to remove
Combo1 <- left_join(Combo, Sample, by = "Samp_mass")
Combo2 <- left_join(Combo1, Blank, by = "Blank_mass")

#Now we can calculate the error between the masses in the sample and those in the blank

Combo3 <- Combo2 %>% mutate(err_ppm = ((Samp_mass-Blank_mass)/(Samp_mass))*10^6)

#You can decide what error limit you want to allow, though 3 ppm is typically pretty good.
#First you can separate it into a dataframe of sample masses that are unique relative to the
#blank and those that are common in both.

Unique <- Combo3 %>% filter(abs(err_ppm) > 3)
Common <- Combo3 %>% filter(abs(err_ppm) <= 3)

#From this point you can either go forward with only the ones unique to the sample and leave
#all common ones behind, or you can go ratio or subtraction evaluation described in the previous 
#section to see which of the common ones you want to retain.

#You will also need to clean up the dataframes to only have two (or three in case of LC) columns 
#before feeding the mass list into the MFAssignR functions.

#All of this can be tweaked as desired, but it provides a decent framework for how I typically
#work through background subtraction.