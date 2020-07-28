						*****************************
						*	Norma Arrazola Herrera	*
						*			Thesis	 		*
						*	 Mediation Analysis		*
						*		  Human CEA	 		*
						*		  Versión 1			*
						*		 5 July 2020	 	*
						*****************************

//Following tutorial from https://stats.idre.ucla.edu/stata/faq/how-can-i-analyze-multiple-mediators-in-stata/
						
clear
set more off
use "InheritedDigitalDivide_sample.dta",clear
numlabel _all, add
count	


************		MEDIATION ANALYSIS FOR HUMAN CEA			************
*1) First we obtain the total effect c with a simple regression of eduparent on HCEA (equation 1 in Appendix G): 
*We start with a simple regression to obtain the total effect (c)
reg HCEA_sd edu_parent_sd sexox_sd edadx_sd residence_type_sd

*2) Then using sureg we obtain a’s, b’s and c’: 
// this will calculate the 4 regression equations (3 times equation 3 (for each mediator) of Appendix G and the equation 2 of Appendix G)
sureg (edu_userSD_sd edu_parent_sd sexox_sd edadx_sd residence_type_sd) (skills_mri_sd edu_parent_sd sexox_sd edadx_sd residence_type_sd) (ict_requipment_sd edu_parent_sd sexox_sd edadx_sd residence_type_sd) (HCEA_sd edu_parent_sd edu_userSD_sd skills_mri_sd ict_requipment_sd sexox_sd edadx_sd residence_type_sd)

*3) Next we calculate the ab for each mediator and for the sum of the three mediators using lncom postestimation function for each mediator
// sureg does not retun the ab term, so we use nl com.
// First the indirect effect through:
*3.1) users education:
nlcom [edu_userSD_sd]_b[edu_parent_sd] * [HCEA_sd]_b[edu_userSD_sd]

*3.2) skills:
nlcom [skills_mri_sd]_b[edu_parent_sd] * [HCEA_sd]_b[skills_mri_sd]

*3.3) ict_equipment:
nlcom [ict_requipment_sd]_b[edu_parent_sd] * [HCEA_sd]_b[ict_requipment_sd]

*3.4 Total indirect effect
// we can add up our indirect effects to get the total indirect effect:

nlcom [edu_userSD_sd]_b[edu_parent_sd] * [HCEA_sd]_b[edu_userSD_sd] + [ict_requipment_sd]_b[edu_parent_sd] * [HCEA_sd]_b[ict_requipment_sd] + [skills_mri_sd]_b[edu_parent_sd] * [HCEA_sd]_b[skills_mri_sd]


*4) Next the bootstrapped CI’s biased corrected: 
// If we want to use bootstrapping we first have to write our own function
// that stores the required ab statistic after running sureg, so that the 
// bootstrap function can do bootstrapping on these ab as obtained with 1000 or
// more bootstrap functions


capture program drop multimed

program multimed, rclass
	sureg (edu_userSD_sd edu_parent_sd sexox_sd edadx_sd residence_type_sd) (skills_mri_sd edu_parent_sd sexox_sd edadx_sd residence_type_sd) (ict_requipment_sd edu_parent_sd sexox_sd edadx_sd residence_type_sd) (HCEA_sd edu_parent_sd edu_userSD_sd skills_mri_sd ict_requipment_sd sexox_sd edadx_sd residence_type_sd)

  return scalar indedu   = [edu_userSD_sd]_b[edu_parent_sd] * [HCEA_sd]_b[edu_userSD_sd]
  
  return scalar indequ   = [ict_requipment_sd]_b[edu_parent_sd] * [HCEA_sd]_b[ict_requipment_sd]
  
  return scalar indski   = [skills_mri_sd]_b[edu_parent_sd] * [HCEA_sd]_b[skills_mri_sd]
  
  return scalar indtotal = [edu_userSD_sd]_b[edu_parent_sd] * [HCEA_sd]_b[edu_userSD_sd] + [ict_requipment_sd]_b[edu_parent_sd] * [HCEA_sd]_b[ict_requipment_sd] + [skills_mri_sd]_b[edu_parent_sd] * [HCEA_sd]_b[skills_mri_sd]

end

// next perform bootstrapping on the newly created function
bootstrap r(indedu) r(indequ) r(indski) r(indtotal), reps(1000): multimed
estat bootstrap, bc

*5) Finally, we bootstrap the differences between the ab of the various mediators in a pairwise fashion to obtain bias corrected CI for each of these differences: 
capture program drop pwcmed

program pwcmed, rclass
  	sureg (edu_userSD_sd edu_parent_sd sexox_sd edadx_sd residence_type_sd) (skills_mri_sd edu_parent_sd sexox_sd edadx_sd residence_type_sd) (ict_requipment_sd edu_parent_sd sexox_sd edadx_sd residence_type_sd) (HCEA_sd edu_parent_sd edu_userSD_sd skills_mri_sd ict_requipment_sd sexox_sd edadx_sd residence_type_sd)

  
  return scalar eduminequ   = [edu_userSD_sd]_b[edu_parent_sd] * [HCEA_sd]_b[edu_userSD_sd] - [ict_requipment_sd]_b[edu_parent_sd] * [HCEA_sd]_b[ict_requipment_sd]
  
  return scalar eduminski   = [edu_userSD_sd]_b[edu_parent_sd] * [HCEA_sd]_b[edu_userSD_sd] - [skills_mri_sd]_b[edu_parent_sd] * [HCEA_sd]_b[skills_mri_sd]
  
  return scalar skiminequ   = [skills_mri_sd]_b[edu_parent_sd] * [HCEA_sd]_b[skills_mri_sd] - [ict_requipment_sd]_b[edu_parent_sd] * [HCEA_sd]_b[ict_requipment_sd]
  
end

bootstrap r(eduminequ) r(eduminski) r(skiminequ), reps(1000): pwcmed
estat bootstrap, bc
