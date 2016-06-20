# Noncoding Mutation Recurrence Prediction

This project is about exploration of regression models for noncoding mutation recurrence in cancer, aiming at finding a most accurate regression model for noncoding mutation recurrence prediction so that based on the predicted recurrences, we can rank the noncoding mutations by their recurrences from high to low, and then choose the top N mutations as cancer driver mutation candidates and then identify through biological experiment validations. 


<b>Author:</b> Tanjin Xu, tjx711@gmail.com <br>
<b>Respondent:</b> Dr. Stephen A. Ramsey, stephen.ramsey@oregonstate.edu <br>
<b>Org:</b> Dr. Ramsey Lab, Oregon State University, http://lab.saramsey.org/

<img src="overall-proc.png">

<p> 
<b>Mainly it includes three steps in the overall procedure:</b>
<ol>
<li> Noncoding mutation annotation, feature extraction, and recurrence calculation </li>
<li> Regression analysis with three main categories of nonlinear models: 
   <ul>
     <li> Generalized linear model including Poisson and Negative Binomial (written in R) </li>
     <li> Ensemble of decision trees including Random Forest and Boosting (written in Python)</li>
     <li> Deep neural network (written in Python)</li>
   </ul>
</li>
<li> Result analysis and plotting </li>
</ol>
</p>

<p>
<b> Project directory structure: </b>
</p>

<p>
<b>Code examples explanation:</b> 
<ul>
<li> <a href="https://github.com/tj711/noncoding-recurrence/blob/master/source/noncoding_extract_features.py"> noncoding_extract_features.py </a> : the main entry to extract features.  </li>
<li> <a href="https://github.com/tj711/noncoding-recurrence/blob/master/source/noncoding_mut_stats.py"> noncoding_mut_stats.py </a> : a program to calculate the mutation recurrences within 101 bps window.</li>
<li> <a href="https://github.com/tj711/noncoding-recurrence/blob/master/source/run_features.sh"> run_features.sh </a> : the bash script to call feature extracting with different inputs of raw features in BED format. </li>
<li> <a href="https://github.com/tj711/noncoding-recurrence/blob/master/source/run_stats.sh"> run_stats.sh </a> : a bash scripts to call the calculation procedure of noncoding mutation recurrence. </li>
<li> <a href="https://github.com/tj711/noncoding-recurrence/blob/master/source/join_features.sh"> join_features.sh </a> : the bash script to join all the extracted features to a big matrix. </li>
<li> <a href="https://github.com/tj711/noncoding-recurrence/blob/master/source/noncoding_mut_anno.py"> noncoding_mut_anno.py </a> : a program to annotate the noncoding mutations out of all the input mutations downloaded from the COSMIC database. </li>
<li> <a href="https://github.com/tj711/noncoding-recurrence/blob/master/source/noncoding_mut_regression.py"> noncoding_mut_regression.py </a> : the main file for exploration of conventional machine learning regression models including randomforest, adaboost, gradientBoosting.</li>
<li> <a href="https://github.com/tj711/noncoding-recurrence/blob/master/source/noncoding_glm_regr.v4.R"> noncoding_glm_regr.v4.R </a> : the main file for exploration of generalized linear model including Poisson and Negative Binomial.</li>
<li> <a href="https://github.com/tj711/noncoding-recurrence/blob/master/source/noncoding_sae_baseline.py"> noncoding_sae_baseline.py </a> : a simple artifical neural network with one hidden layer regression.</li>
<li> <a href="https://github.com/tj711/noncoding-recurrence/blob/master/source/noncoding_sae.py"> noncoding_sae.py </a> : stacked auto encoder regression.</li>
</ul>
</p>

<p>
<b>Remaining work:</b>
<ul>
<li> <b>Add more features.</b> Some considerations: 1) introduce more specific Transcription Factors; 2) the distance to the oncogenes/tumor-suppressor genes; 3) the DNA shape/histone info; 4) chromosome length; 5) etc..</li>
<li> <b>Left work for regression analysis.</b> 1) fit zero-truncated negative binomial regression into our current noncoding mutation data; 2) need special handling with the super outliers with rare extreme high frequency; 3) check the predicted recurrence of the Poisson/NB model, i.e., check the top 100 mutations with highest predicted recurrence in comparison to the true frequency value; 4) maybe a ranking regressoin problem; 5) need carefully re-lookinto the features data. </li>
<li> <b>Binary classfication problem.</b> The challenge is to construct balanced training set which is consist of the positive cases (true noncoding mutation drivers) and the negative cases (validated passenger mutations). </li>
<li> <b>Focus on specific categories of cancer types.</b> Considering overall cancer types together may be more difficult due to high heterogeneity of different tumor samples in different cancer types. </li>
</ul>
</p>
