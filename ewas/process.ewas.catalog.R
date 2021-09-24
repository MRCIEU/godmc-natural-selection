library(data.table)
library(tidyverse)
rmeta<-fread("ewascatalog-studies.txt.gz")
r<-fread("ewascatalog-results.txt.gz")
r<-r[r$P<1e-7,]
length(unique(r$StudyID))
#[1] 692 #1998
nrow(data.frame(table(rmeta$Exposure)))
#[1] 647 #1835
studyid<-unique(r$StudyID)
rmeta<-rmeta[which(rmeta$StudyID%in%studyid),]
length(unique(rmeta$StudyID)) #1998
length(unique(r$StudyID)) #1998

g1<-grep("whole blood",rmeta$Tissue,ignore.case =T)
table(rmeta[g1,"Tissue"])
rmeta<-rmeta[g1,]
studyid<-unique(rmeta$StudyID)
r<-r[which(r$StudyID%in%studyid),]

length(unique(rmeta$StudyID))
#1804
length(unique(r$StudyID))
#[1] 1804

g1<-grep("european",rmeta$Ancestry,ignore.case =T)
g2<-grep("Unclear",rmeta$Ancestry,ignore.case =T)
g1<-unique(c(g1,g2))
table(rmeta[g1,"Ancestry"])
rmeta<-rmeta[g1,]
studyid<-unique(rmeta$StudyID)
r<-r[which(r$StudyID%in%studyid),]
length(unique(r$StudyID))
# [1] 1754
length(unique(rmeta$StudyID))
#[1] 1754

rmeta<-rmeta[rmeta$N>1000,]
studyid<-unique(rmeta$StudyID)
r<-r[which(r$StudyID%in%studyid),]
studyid<-unique(r$StudyID)
rmeta<-rmeta[which(rmeta$StudyID%in%studyid),]
length(unique(rmeta$StudyID))
#1462
length(unique(r$StudyID))
#1462
nrow(data.frame(table(rmeta$Exposure)))
#1203
nrow(data.frame(table(rmeta$Outcome)))
#114

nrow(r)
nrow(rmeta)
m<-match(r$StudyID,rmeta$StudyID)
rmeta2<-data.frame(rmeta[m,])
nrow(rmeta2)
r<-data.frame(r)

r<-data.frame(r,Outcome=rmeta2[,"Outcome"],Exposure=rmeta2[,"Exposure"])
temp<-data.frame(Outcome=as.character(r$Outcome),Exposure=as.character(r$Exposure))
temp$Outcome<-as.character(gsub("DNA methylation",NA,temp$Outcome))
temp$Exposure<-as.character(gsub("DNA methylation",NA,temp$Exposure))
temp[is.na(temp)] = ''
temp<-unite(temp, Phenotype, 1:2, sep='')
r$Phenotype<-temp$Phenotype

mondo<-fread("~/ARIES/methylation/ICC_2018/data/mondo_efo_mappings.tsv",he=F)
mondo$V2<-gsub("http://www.ebi.ac.uk/efo/","",mondo$V2)
mondo$V1<-gsub("http://purl.obolibrary.org/obo/","",mondo$V1)

out<-data.frame(mondo_id=1,mondo_name=2)
for (i in 1:nrow(rmeta)){
m<-match(rmeta[i,"EFO"],mondo$V2)
out[i,1]<-paste0(mondo[m,V1])
out[i,2]<-paste0(mondo[m,V3])
}
rmeta$mondo_id<-out[,1]
rmeta$mondo_name<-out[,2]

#ariesid<-rmeta[which(rmeta$Consortium=="ARIES"),]
#ariesid$StudyID
#aries<-r[which(r$StudyID%in%ariesid$StudyID),]
#150
#length(unique(aries$StudyID))
#99
r$subcategory<-NA

g<-grep("cholesterol",r$StudyID)
r$subcategory[g]<-"Lipid"
g<-grep("hdl",r$StudyID)
r$subcategory[g]<-"Lipid"
g<-grep("ldl",r$StudyID)
r$subcategory[g]<-"Lipid"
g<-grep("triglycerides",r$StudyID)
r$subcategory[g]<-"Lipid"
g<-grep("idl",r$StudyID)
r$subcategory[g]<-"Lipid"
g<-grep("sphingomyelins",r$StudyID)
r$subcategory[g]<-"Lipid"
g<-grep("sphingomyelins",r$StudyID)
r$subcategory[g]<-"Lipid"
g<-grep("highdensity_lipoprotein",r$StudyID)
r$subcategory[g]<-"Lipid"
g<-grep("lowdensity_lipoprotein",r$StudyID)
r$subcategory[g]<-"Lipid"
g<-grep("lipemia",r$StudyID)
r$subcategory[g]<-"Lipid"
g<-grep("hypertriglyceridemic_waist",r$StudyID)
r$subcategory[g]<-"Lipid"
g<-grep("statin_use",r$StudyID)
r$subcategory[g]<-"Lipid"
g<-grep("32237968_Jones-GT_lipoprotein_a",r$StudyID)
r$subcategory[g]<-"Lipid"
g<-grep("Lp(a)",r$Phenotype,fixed=T)
r$subcategory[g]<-"Lipid"

#g<-grep("fatty_acid",r$StudyID)
#r$subcategory[g]<-"Fatty acid"
#g<-grep("pufa",r$StudyID)
#r$subcategory[g]<-"Fatty acid"
#g<-grep("arachidonoylglycerophosphoethanolamine",r$StudyID)
#r$subcategory[g]<-"Metabolite"
#g<-grep("palmitoylglycerophosphocholine",r$StudyID)
#r$subcategory[g]<-"Metabolite"
#g<-grep("butyrylcarnitine",r$StudyID)
#r$subcategory[g]<-"Metabolite"
#g<-grep("indoxyl_sulfate",r$StudyID)
#r$subcategory[g]<-"Metabolite"
#g<-grep("palmitoleate",r$StudyID)
#r$subcategory[g]<-"Metabolite"
#g<-grep("24014485_Petersen",r$StudyID)
#r$subcategory[g]<-"Metabolite"

#g<-grep("24014485_Petersen",r$StudyID)
#r$subcategory[g]<-"Metabolite"
#g<-grep("27839851_ppdde",r$StudyID)
#r$subcategory[g]<-"Metabolite"


g<-grep("29740313_Mendelson-MM_soluble_tumor_necrosis_factor_receptor_2",r$StudyID)
r$subcategory[g]<-"Cytokine"
g<-grep("tumour_necrosis_factor",r$StudyID)
r$subcategory[g]<-"Cytokine"

g<-grep("vascular_endothelial_growth_factor_a",r$StudyID)
r$subcategory[g]<-"VEGFA"


g<-grep("atrial_fibrillation",r$StudyID)
r$subcategory[g]<-"Cardiovascular"
g<-grep("myocardial_infarction",r$StudyID)
r$subcategory[g]<-"Cardiovascular"
g<-grep("myocardial_infarction",r$StudyID)
r$subcategory[g]<-"Cardiovascular"
g<-grep("coronary_heart_disease",r$StudyID)
r$subcategory[g]<-"Cardiovascular"
g<-grep("cardiovascular_risk",r$StudyID)
r$subcategory[g]<-"Cardiovascular"
g<-grep("stroke",r$StudyID)
r$subcategory[g]<-"Cardiovascular"

g<-grep("hypertensive",r$StudyID)
r$subcategory[g]<-"Hypertension"
g<-grep("blood_pressure",r$StudyID)
r$subcategory[g]<-"Hypertension"

g<-grep("gingival",r$StudyID)
r$subcategory[g]<-"Teeth"
g<-grep("tooth",r$StudyID)
r$subcategory[g]<-"Teeth"

g<-grep("skeletal_muscle_fiber_type",r$StudyID)
r$subcategory[g]<-"Skeletal Muscle"


#g<-grep("pyruvate",r$StudyID)
#r$subcategory[g]<-"Carbohydrate"
#g<-grep("glycoprotein_acetyls",r$StudyID)
#r$subcategory[g]<-"Protein"

g<-grep("c-reactive_protein",r$StudyID)
r$subcategory[g]<-"Immune system"
g<-grep("creactive_protein",r$StudyID)
r$subcategory[g]<-"Immune system"
g<-grep("inflammatory_biomarker",r$StudyID)
r$subcategory[g]<-"Immune system"


g<-grep("waist_circumference",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("fat_mass",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("region_fat",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("tissue_fat",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("body_mass_index",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("underweight",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("overweight",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("hip_circumference",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("obesity",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("waist_circumfrence",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("subcutaneous_fat_thickness",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("bmi_discovery",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("Battram-T_weight",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("skinfold_thickness",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("arm_circumference",r$StudyID)
r$subcategory[g]<-"Anthropometric"
g<-grep("subscapular_skinfold_thickness",r$StudyID)
r$subcategory[g]<-"Anthropometric"
#g<-grep("mass",r$StudyID)
#r$subcategory[g]<-"Anthropometric"


g<-grep("adiponectin",r$StudyID)
r$subcategory[g]<-"Adiponectin/Resistin"
g<-grep("plasma_resistin",r$StudyID)
r$subcategory[g]<-"Adiponectin/Resistin"


g<-grep("insulin",r$StudyID)
r$subcategory[g]<-"Glycemic"
g<-grep("diabetes",r$StudyID)
r$subcategory[g]<-"Glycemic"
g<-grep("glucose",r$StudyID)
r$subcategory[g]<-"Glycemic"
g<-grep("homa",r$StudyID)
r$subcategory[g]<-"Glycemic"


g<-grep("iron",r$StudyID)
r$subcategory[g]<-"Metals"
g<-grep("nickel",r$StudyID)
r$subcategory[g]<-"Metals"
g<-grep("arsenic",r$StudyID)
r$subcategory[g]<-"Metals"
g<-grep("selenium",r$StudyID)
r$subcategory[g]<-"Metals"
g<-grep("cadmium",r$StudyID)
r$subcategory[g]<-"Metals"
g<-grep("vanadium",r$StudyID)
r$subcategory[g]<-"Metals"
g<-grep("copper",r$StudyID)
r$subcategory[g]<-"Metals"
g<-grep("lead",r$StudyID)
r$subcategory[g]<-"Metals"

g<-grep("breastfeeding",r$StudyID)
r$subcategory[g]<-"Baby Feeding"
g<-grep("bottle_feeding",r$StudyID)
r$subcategory[g]<-"Baby Feeding"

g<-grep("kidney",r$StudyID)
r$subcategory[g]<-"Kidney"

g<-grep("liver",r$StudyID)
r$subcategory[g]<-"Liver"
g<-grep("gammaglutamyl",r$StudyID)
r$subcategory[g]<-"Liver"
g<-grep("hepatic_fat",r$StudyID)
r$subcategory[g]<-"Liver"
g<-grep("aspartate_aminotransferase",r$StudyID)
r$subcategory[g]<-"Liver"
g<-grep("alanine_aminotransferase",r$StudyID)
r$subcategory[g]<-"Liver"



g<-grep("lung_function",r$StudyID)
r$subcategory[g]<-"Lung function"
g<-grep("fev1fvc",r$StudyID)
r$subcategory[g]<-"Lung function"
g<-grep("fev1",r$StudyID)
r$subcategory[g]<-"Lung function"
g<-grep("fvc",r$StudyID)
r$subcategory[g]<-"Lung function"
g<-grep("diffusing_capacity_of_the_lung_for_carbon_monoxide",r$StudyID)
r$subcategory[g]<-"Lung function"
g<-grep("copd",r$StudyID)
r$subcategory[g]<-"Lung function"


g<-grep("adversity",r$StudyID)
r$subcategory[g]<-"Adversity"
g<-grep("child_abuse",r$StudyID)
r$subcategory[g]<-"Adversity"
g<-grep("perceived_discrimination",r$StudyID)
r$subcategory[g]<-"Adversity"
g<-grep("domestic_violence",r$StudyID)
r$subcategory[g]<-"Adversity"
g<-grep("sexual_victimization",r$StudyID)
r$subcategory[g]<-"Adversity"
g<-grep("bullying",r$StudyID)
r$subcategory[g]<-"Adversity"
g<-grep("sexual_abuse",r$StudyID)
r$subcategory[g]<-"Adversity"


g<-grep("bone",r$StudyID)
r$subcategory[g]<-"Bone mineral density"

#g<-grep("acetate",r$StudyID)
#r$subcategory[g]<-"AminoAcid"
#g<-grep("homocysteine",r$StudyID)
#r$subcategory[g]<-"AminoAcid"


g<-grep("28556291_prenatal_phthalate_exposure_26_weeks_gestation",r$StudyID)
r$subcategory[g]<-"Chemicals"
g<-grep("aflatoxin_b1",r$StudyID)
r$subcategory[g]<-"Chemicals"

g<-grep("folate",r$StudyID)
r$subcategory[g]<-"Vitamins"
g<-grep("vitamin_d",r$StudyID)
r$subcategory[g]<-"Vitamins"

g<-grep("sex_hormone",r$StudyID)
r$subcategory[g]<-"Sex hormone"
g<-grep("sex_hormone",r$StudyID)
r$subcategory[g]<-"Sex hormone"
g<-grep("esterogen",r$StudyID)
r$subcategory[g]<-"Sex hormone"

g<-grep("hiv_infection",r$StudyID)
r$subcategory[g]<-"Infectious Disease"
g<-grep("human_immunodeficiency_virus",r$StudyID)
r$subcategory[g]<-"Infectious Disease"
g<-grep("acquired_hiv",r$StudyID)
r$subcategory[g]<-"Infectious Disease"
g<-grep("hepatitis_c",r$StudyID)
r$subcategory[g]<-"Infectious Disease"

g<-grep("carcinoma",r$StudyID)
r$subcategory[g]<-"Cancer"
g<-grep("cancer",r$StudyID)
r$subcategory[g]<-"Cancer"
g<-grep("adenomas",r$StudyID)
r$subcategory[g]<-"Cancer"
g<-grep("melanoma",r$StudyID)
r$subcategory[g]<-"Cancer"
g<-grep("squamous_cell",r$StudyID)
r$subcategory[g]<-"Cancer"

g<-grep("age_and_sex",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("_age_t_cells",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("_age_monocytes",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("age_4_vs_age_0",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("28811542_age",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("25888029_ageing",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("24163245_age",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("28693600_age",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("age_ewas_catalog",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("Tajuddin-SM_age",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("31892350_McCartney-DL_age_replication",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("31892350_McCartney-DL_age_discovery",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("age_model",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("age_fetal",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("age_adult",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("Wang-Y_aging",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("Jansen-R_age",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("leukocyte_telomere_length",r$StudyID)
r$subcategory[g]<-"Aging"
g<-grep("mortality",r$StudyID)
r$subcategory[g]<-"Aging"


g<-grep("migration",r$StudyID)
r$subcategory[g]<-"Geography"
g<-grep("ethnicity",r$StudyID)
r$subcategory[g]<-"Geography"
g<-grep("altitude",r$StudyID)
r$subcategory[g]<-"Geography"
g<-grep("papuan_ancestry",r$StudyID)
r$subcategory[g]<-"Geography"


g<-grep("air_pollution",r$StudyID)
r$subcategory[g]<-"Air pollution"
g<-grep("particulate_matter",r$StudyID)
r$subcategory[g]<-"Air pollution"
g<-grep("no2_exposure",r$StudyID)
r$subcategory[g]<-"Air pollution"
g<-grep("pm2_5",r$StudyID)
r$subcategory[g]<-"Air pollution"
g<-grep("Chi-GC_nox",r$StudyID)
r$subcategory[g]<-"Air pollution"
g<-grep("nitrogen_dioxide",r$StudyID)
r$subcategory[g]<-"Air pollution"

g<-grep("gestational_age",r$StudyID)
r$subcategory[g]<-"Reproductive"
g<-grep("Kashima-K_birth_weight",r$StudyID)
r$subcategory[g]<-"Reproductive"
g<-grep("Antoun-E_birthweight",r$StudyID)
r$subcategory[g]<-"Reproductive"
g<-grep("Hannon-E_birthweight",r$StudyID)
r$subcategory[g]<-"Reproductive"
g<-grep("Kupers-L_birthweight",r$StudyID)
r$subcategory[g]<-"Reproductive"
g<-grep("birth_weight_twinsuk",r$StudyID)
r$subcategory[g]<-"Reproductive"
g<-grep("28264723_birth_weight",r$StudyID)
r$subcategory[g]<-"Reproductive"
g<-grep("24561991_birth_weight",r$StudyID)
r$subcategory[g]<-"Reproductive"
g<-grep("31600381_Kresovich-JK_parity",r$StudyID)
r$subcategory[g]<-"Reproductive"

g<-grep("rapid_weight_gain",r$StudyID)
r$subcategory[g]<-"Reproductive"

g<-grep("28165855_early_spontaneous_preterm_birth_discovery",r$StudyID)
r$subcategory[g]<-"Reproductive"
g<-grep("pubertal_timing",r$StudyID)
r$subcategory[g]<-"Reproductive"
g<-grep("preeclampsia",r$StudyID)
r$subcategory[g]<-"Reproductive"
g<-grep("fetal_intolerance_of_labor",r$StudyID)
r$subcategory[g]<-"Reproductive"
g<-grep("age_at_menarche",r$StudyID)
r$subcategory[g]<-"Reproductive"
g<-grep("season_of_birth",r$StudyID)
r$subcategory[g]<-"Reproductive"


g<-grep("sex_kora",r$StudyID)
r$subcategory[g]<-"Sex"
g<-grep("sex_discovery",r$StudyID)
r$subcategory[g]<-"Sex"
g<-grep("sex_replication",r$StudyID)
r$subcategory[g]<-"Sex"
g<-grep("26500701_sex_metaanalysis",r$StudyID)
r$subcategory[g]<-"Sex"
g<-grep("26553366_sex",r$StudyID)
r$subcategory[g]<-"Sex"
g<-grep("25650246_sex",r$StudyID)
r$subcategory[g]<-"Sex"
g<-grep("sex_ewas_catalog",r$StudyID)
r$subcategory[g]<-"Sex"

g<-grep("age-by-sex_interaction",r$StudyID)
r$subcategory[g]<-"Age by sex interaction"

g<-grep("rheumatoid_arthritis",r$StudyID)
r$subcategory[g]<-"Autoimmune / inflammatory"
g<-grep("sjogrens_syndrome",r$StudyID)
r$subcategory[g]<-"Autoimmune / inflammatory"
g<-grep("lupus",r$StudyID)
r$subcategory[g]<-"Autoimmune / inflammatory"
g<-grep("inflammatory_bowel_disease",r$StudyID)
r$subcategory[g]<-"Autoimmune / inflammatory"
g<-grep("crohn",r$StudyID)
r$subcategory[g]<-"Autoimmune / inflammatory"
g<-grep("ulcerative_colitis",r$StudyID)
r$subcategory[g]<-"Autoimmune / inflammatory"
g<-grep("graves",r$StudyID)
r$subcategory[g]<-"Autoimmune / inflammatory"
g<-grep("multiple_sclerosis",r$StudyID)
r$subcategory[g]<-"Autoimmune / inflammatory"


g<-grep("psoriasis",r$StudyID)
r$subcategory[g]<-"Skin disease"
g<-grep("sunlight",r$StudyID)
r$subcategory[g]<-"Skin disease"
g<-grep("total_body_naevus_count",r$StudyID)
r$subcategory[g]<-"Skin disease"
g<-grep("nevus_count",r$StudyID)
r$subcategory[g]<-"Skin disease"

g<-grep("triiodothyronine",r$StudyID)
r$subcategory[g]<-"Thyroid function"
g<-grep("thyrotropin",r$StudyID)
r$subcategory[g]<-"Thyroid function"

g<-grep("atopy",r$StudyID)
r$subcategory[g]<-"Allergy"
g<-grep("asthma",r$StudyID)
r$subcategory[g]<-"Respiratory"
g<-grep("allergy",r$StudyID)
r$subcategory[g]<-"Allergy"

g<-grep("total_serum_ige",r$StudyID)
r$subcategory[g]<-"Immune cell subset"
g<-grep("growth_differentiation_factor_15",r$StudyID)
r$subcategory[g]<-"Immune cell subset"

g<-grep("28535255_tea_consumption_metaanalysis",r$StudyID)
r$subcategory[g]<-"Diet"
g<-grep("fruit_consumption",r$StudyID)
r$subcategory[g]<-"Diet"
g<-grep("mediterranean_diet",r$StudyID)
r$subcategory[g]<-"Diet"
g<-grep("healthy_eating",r$StudyID)
r$subcategory[g]<-"Diet"
g<-grep("low-energy_diet",r$StudyID)
r$subcategory[g]<-"Diet"

g<-grep("cleft",r$StudyID)
r$subcategory[g]<-"Birth defects"

g<-grep("schizophrenia",r$StudyID)
r$subcategory[g]<-"Psychiatric disorders"
g<-grep("depression",r$StudyID)
r$subcategory[g]<-"Psychiatric disorders"
g<-grep("depressive",r$StudyID)
r$subcategory[g]<-"Psychiatric disorders"
g<-grep("28595673_earlyonset_conduct_problems",r$StudyID)
r$subcategory[g]<-"Psychiatric disorders"
g<-grep("attention_deficit",r$StudyID)
r$subcategory[g]<-"Psychiatric disorders"
g<-grep("dementia",r$StudyID)
r$subcategory[g]<-"Psychiatric disorders"

g<-grep("amyloid",r$StudyID)
r$subcategory[g]<-"Neurological disorders"
g<-grep("alzheimer",r$StudyID)
r$subcategory[g]<-"Neurological disorders"
g<-grep("wellbeing",r$StudyID)
r$subcategory[g]<-"Neurological disorders"
g<-grep("parkinson",r$StudyID)
r$subcategory[g]<-"Neurological disorders"
g<-grep("palsy",r$StudyID)
r$subcategory[g]<-"Neurological disorders"
g<-grep("neurodegenerative_disorders",r$StudyID)
r$subcategory[g]<-"Neurological disorders"
g<-grep("seizure",r$StudyID)
r$subcategory[g]<-"Neurological disorders"



g<-grep("sleepiness",r$StudyID)
r$subcategory[g]<-"Psychiatric disorders"
g<-grep("neurobehavioural",r$StudyID)
r$subcategory[g]<-"Psychiatric disorders"
g<-grep("antidepressant",r$StudyID)
r$subcategory[g]<-"Psychiatric disorders"
g<-grep("anxiety",r$StudyID)
r$subcategory[g]<-"Psychiatric disorders"
g<-grep("psychosis",r$StudyID)
r$subcategory[g]<-"Psychiatric disorders"
g<-grep("aggressive_behaviour",r$StudyID)
r$subcategory[g]<-"Psychiatric disorders"
g<-grep("borderline",r$StudyID)
r$subcategory[g]<-"Psychiatric disorders"


g<-grep("downs_syndrome",r$StudyID)
r$subcategory[g]<-"Neurodevelopmental disorders"

g<-grep("hba1c",r$StudyID)
r$subcategory[g]<-"Blood"
g<-grep("coagulation_factor_viii",r$StudyID)
r$subcategory[g]<-"Blood"

g<-grep("glutathione",r$StudyID)
r$subcategory[g]<-"Oxidative Stress"
g<-grep("conjugated_dienes",r$StudyID)
r$subcategory[g]<-"Oxidative Stress"
g<-grep("homocysteine",r$StudyID)
r$subcategory[g]<-"Oxidative Stress"

g<-grep("smoking",r$StudyID)
r$subcategory[g]<-"Behavioural"
g<-grep("alcohol",r$StudyID)
r$subcategory[g]<-"Behavioural"
g<-grep("cotinine",r$StudyID)
r$subcategory[g]<-"Behavioural"
g<-grep("physical_activity",r$StudyID)
r$subcategory[g]<-"Behavioural"
g<-grep("substance_use",r$StudyID)
r$subcategory[g]<-"Behavioural"
g<-grep("smoke_exposure",r$StudyID)
r$subcategory[g]<-"Behavioural"



g<-grep("socioeconomic_position",r$StudyID)
r$subcategory[g]<-"Socioeconomic Position"
g<-grep("educational_attainment",r$StudyID)
r$subcategory[g]<-"Education"
g<-grep("education",r$StudyID)
r$subcategory[g]<-"Education"
g<-grep("cognitive_abilities",r$StudyID)
r$subcategory[g]<-"Education"

g<-grep("31062658_Santos",r$StudyID)
r$subcategory[g]<-"Social"

g<-grep("thalamus_volume",r$StudyID)
r$subcategory[g]<-"Brain trait"
g<-grep("fetal_brain",r$StudyID)
r$subcategory[g]<-"Brain trait"

g<-grep("30602389_Battram-T_tissue_ewas_catalog",r$StudyID)
r$subcategory[g]<-"Tissue"

g<-grep("31282290_Huan-T_mir",r$StudyID)
r$subcategory[g]<-"mi-RNA"

g<-grep("fatty_acid",r$StudyID)
r$subcategory[g]<-"Fatty acid"
g<-grep("pufa",r$StudyID)
r$subcategory[g]<-"Fatty acid"

g<-grep("24014485_Petersen",r$StudyID)
r$subcategory[g]<-"Metabolite"

g<-grep("33413638_Gomez",r$StudyID)
r$subcategory[g]<-"Metabolite"

g<-grep("29055820_Wahl",r$StudyID)
r$subcategory[g]<-"IgG glycosylation"


g<-grep("33752734_Yao-C",r$StudyID)
r$subcategory[g]<-"Gene expression"

g<-grep("33550919_Huang-R",r$StudyID)

df.na<-data.frame(r[which(is.na(r$subcategory)),c("StudyID","subcategory")])
unique(df.na$StudyID)

length(r$CpG)
#[1] 546998
length(unique(r$CpG))
#[1] 312820

save(r,rmeta,file="ewas_category.Robj")

table(r$subcategory)

#             Adversity Age by sex interaction                  Aging 
#                    59                     22                 470614 
#         Air pollution         Anthropometric            Behavioural 
#                  7262                   1473                  14751 
#  Bone mineral density            Brain trait                 Cancer 
#                     3                      2                   1872 
#        Cardiovascular               Cytokine                   Diet 
#                    39                    175                     25 
#             Education        Gene expression               Glycemic 
#                   166                   3418                    127 
#          Hypertension      IgG glycosylation     Immune cell subset 
#                   160                      6                     13 
#         Immune system                 Kidney                  Lipid 
#                  217                     99                    648 
#                 Liver          Lung function             Metabolite 
#                    51                    125                   2673 
#                Metals                 mi-RNA Neurological disorders 
#                     4                    123                     13 
#      Oxidative Stress  Psychiatric disorders           Reproductive 
#                   106                   2301                   4328 
#           Respiratory                    Sex       Thyroid function 
#                   41                  36052                      7 
#             Vitamins 
#                   2 

#w<-which(r$subcategory%in%c("Metabolite","Gene expression","mi-RNA","Tissue","Brain trait"))
#r<-r[-w,]
nrow(r) #546998
w<-which(r$subcategory%in%c("Adversity","Age by sex interaction","Aging","Air pollution","Behavioural","Diet","Gene expression","Metabolite","IgG glycosylation","Metals","mi-RNA","Oxidative Stress","Sex","Vitamins"))
r2<-r[-w,]
nrow(r2)
#11881

g<-grep("Cancer treatment",r2$Phenotype)
r2<-r2[-g,]
g<-grep("APOE e2 vs e4",r2$Phenotype)
r2<-r2[-g,]
g<-grep("Dichlorodiphenyldichloroethylene",r2$Phenotype)
r2<-r2[-g,]
g<-grep("Statin use",r2$Phenotype)
r2<-r2[-g,]

nrow(r2)
#10487

r2$Phenotype<-gsub("Incident atrial fibrillation","Atrial fibrillation",r2$Phenotype)
r2$Phenotype<-gsub("Prevalent atrial fibrillation","Atrial fibrillation",r2$Phenotype)
r2$Phenotype<-gsub("Lp(a)","Lipoprotein(a) levels",r2$Phenotype,fixed=T)
r2$Phenotype<-gsub("Serum high-density lipoprotein cholesterol","High-density lipoprotein cholesterol",r2$Phenotype)
r2$Phenotype<-gsub("Serum low-density lipoprotein cholesterol","Low-density lipoprotein cholesterol",r2$Phenotype)
r2$Phenotype<-gsub("Serum total cholesterol","Total cholesterol",r2$Phenotype)
r2$Phenotype<-gsub("Serum triglycerides","Triglycerides",r2$Phenotype)
r2$Phenotype<-gsub("Type II diabetes","Type 2 diabetes",r2$Phenotype)
r2$Phenotype<-gsub("Waist circumfrence","Waist circumference",r2$Phenotype)
r2$Phenotype<-gsub("FEV1 / FVC","FEV1/FVC",r2$Phenotype,fixed=T)
r2$Phenotype<-gsub("Cognitive abilities : digit test","Cognitive abilities: digit test",r2$Phenotype,fixed=T)
r2$Phenotype<-gsub("fasting glucose","Fasting glucose",r2$Phenotype)
r2$Phenotype<-gsub("fasting insulin","Fasting insulin",r2$Phenotype)
r2$Phenotype<-gsub("body mass index","Body mass index",r2$Phenotype)


df<-data.frame(table(r2$Phenotype))
nrow(df)
#[1] 56

df2<-data.frame(table(r2$subcategory))
nrow(df2)
# 20

data.frame(df2)
#                     Var1 Freq
#1          Anthropometric 1473
#2    Bone mineral density    3
#3             Brain trait    2
#4                  Cancer  506
#5          Cardiovascular   39
#6                Cytokine  175
#7               Education  166
#8                Glycemic  127
#9            Hypertension  160
#10     Immune cell subset   13
#11          Immune system  217
#12                 Kidney   99
#13                  Lipid  641
#14                  Liver   51
#15          Lung function  125
#16 Neurological disorders   13
#17  Psychiatric disorders 2301
#18           Reproductive 4328
#19            Respiratory   41
#20       Thyroid function    7


outcpg <- r2 %>% 
  filter(!is.na(P)) %>%
  group_by(CpG) %>%
  summarise( 
    nsig=sum(P < 1e-7),
    minP=min(P),
    maxP=max(P),
    pheno_uniq=paste(unique(Phenotype),collapse=";"),
    phenocategory_uniq=paste(unique(subcategory),collapse=";"),
    npheno=length(unique(Phenotype)),
    nphenocategory=length(unique(subcategory))
  )
str(outcpg)
nrow(outcpg)
#7138

outcpg[which(outcpg$nsig>10),]

save(outcpg,file="~/repo/godmc-natural-selection/ewas/outcpg.rdata")

load("/projects/MRC-IEU/research/projects/ieu2/p5/091/working/data/tfbs/tfbsbycpg_summarised.Robj")
names(outcpg)<-c("cpg","n_ewas_sig","ewas_minP","ewas_max_P","ewas_pheno_uniq","ewas_pheno_category_uniq","ewas_npheno","ewas_nphenocategory")
m<-match(tfbs2$cpg,outcpg$cpg)
tfbs2<-data.frame(tfbs2,outcpg[m,])

length(which(is.na(tfbs2$ewas_pheno_category_uniq)))
#[1] 185957
length(which(!is.na(tfbs2$ewas_pheno_category_uniq)))
#[1] 4145

save(tfbs2,file="/projects/MRC-IEU/research/projects/ieu2/p5/091/working/data/tfbs/tfbsbycpg_summarised.Robj")
nrow(tfbs2)

