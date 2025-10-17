"""
calculate and plot the confidence intervals of the individual TE categories in the transcript surroundings
"""


import parse_gff as gff

import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D

from sklearn.preprocessing import PolynomialFeatures
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std


def plot_confidence_intervals(before_filepath:str, after_filepath:str, num_sig_transcripts:int, num_all_transcripts:int, win_len = 200, overlapping_windows = True, all_before_filepath:str = "", all_after_filepath:str = "", all_transcripts:int = 0, filename = "cumulative_repeat_presence_around_transcripts_95_perc_CI", modelstats_filename = "pol_reg_sum.txt", legend = True, plot_white_bg = False):
    """
    make a separate plot for each TE category with the foreground/background and before/after lines in it 
    do the polynomial regressions and see the 95% confidence intervals
    I followed this tutorial: https://ostwalprasad.github.io/machine-learning/Polynomial-Regression-using-statsmodel.html
    """
    plot_transparent_bg = True
    if plot_white_bg:
        plot_transparent_bg = False
    
    before_dict = gff.read_dict_from_file(before_filepath)
    before_dict = { key : [int(v)/num_sig_transcripts*100 for v in value] for key, value in before_dict.items()}
    # show percentages of they are too high
    for key, value in before_dict.items():
        for v in value:
            v_ = int(v)
            perc = v_/num_sig_transcripts*100
            if perc>150:
                print(f"{key} : \t {v_}/{num_sig_transcripts} * 100 = {perc:.2f}%")

    after_dict = gff.read_dict_from_file(after_filepath)
    after_dict = { key : [int(v)/num_sig_transcripts*100 for v in value] for key, value in after_dict.items()}

    all_before_dict = {}
    all_after_dict = {}
    if all_before_filepath !="" and all_after_filepath !="" :
        
        all_before_dict = gff.read_dict_from_file(all_before_filepath)
        all_before_dict = { key : [int(v)/all_transcripts*100 for v in value] for key, value in all_before_dict.items()}
        all_after_dict = gff.read_dict_from_file(all_after_filepath)
        all_after_dict = { key : [int(v)/all_transcripts*100 for v in value] for key, value in all_after_dict.items()}
        
    colors = {
        'Unknown' : "#C1C1C1" , # light grey
        # orange
        'DNA' : "#FF9000" , # Princeton orange
        # green
        'LTR' : "#6E8448" , # reseda green
        'RC' : "#8EA861" , # asparagus 
        # red
        'tRNA' : "#C14953" , # bittersweet shimmer
        'rRNA' : "#D0767E" , # old rose
        'snRNA' : "#7A2A30" , # wine
        # blue 
        'LINE' : "#3476AD" , # UCLA blue
        'SINE': "#72A8D5" , # ruddy blue
        'SINE?': "#72A8D5" , # ruddy blue
        # '' : "#2A618D" , #lapis lazuli
        # dark red-brown
        'Low_complexity' : "#3A3335" , # Jet 
        'Satellite' : "#564D4F" , #Wenge 
        'Simple_repeat' : "#827376" , #Taupe gray
    }

    fs = 25 # set font size

    rep_classes = list(before_dict.keys())
    num_bp = len(before_dict[rep_classes[0]]) 
      
    for rep_class in rep_classes:
        if legend:
            fig, ax = plt.subplots(1, 1, figsize=(23, 10))
        else:
            fig, ax = plt.subplots(1, 1, figsize=(21, 12))

        species = gff.split_at_second_occurrence(before_filepath.split("/")[-1])
        species = species.replace("_", ". ")

        max_percentage = 0

        max_before = max(before_dict[rep_class])
        max_after = max(after_dict[rep_class])
        if max_before>max_percentage:
            max_percentage=max_before
        if max_after>max_percentage:
            max_percentage=max_after

        ## use only every 10th bp for better independence
        # before_dict[rep_class] = before_dict[rep_class][::win_len]
        # all_before_dict[rep_class] = all_before_dict[rep_class][::win_len]
        # after_dict[rep_class] = after_dict[rep_class][::win_len]
        # all_after_dict[rep_class] = all_after_dict[rep_class][::win_len]
        
        ## calculate nonoverlapping window means
        before_df = pd.DataFrame(before_dict[rep_class])
        all_before_df = pd.DataFrame(all_before_dict[rep_class])
        after_df = pd.DataFrame(after_dict[rep_class])
        all_after_df = pd.DataFrame(all_after_dict[rep_class])
        
        if overlapping_windows:
            ## overlapping windows
            legend_title = f"{win_len} bp overlapping windows"
            x_before = range(-num_bp, 0)
            x_after = range(0, num_bp)

            before_df_win = before_df.rolling(window=win_len, min_periods=1).median()
            all_before_df_win = all_before_df.rolling(window=win_len, min_periods=1).median()
            after_df_win = after_df.rolling(window=win_len, min_periods=1).median()
            all_after_df_win = all_after_df.rolling(window=win_len, min_periods=1).median()
        else:
            ## nonoverlapping windows
            legend_title = f"{win_len} bp nonoverlapping windows"
            x_before = range(-num_bp, 0, win_len)
            x_after = range(0, num_bp, win_len)

            before_df_win = before_df.groupby(before_df.index // win_len).median()
            all_before_df_win = all_before_df.groupby(all_before_df.index // win_len).median()
            after_df_win = after_df.groupby(after_df.index // win_len).median()
            all_after_df_win = all_after_df.groupby(all_after_df.index // win_len).median()

        # make to list again
        before_dict[rep_class] = list(before_df_win[0])
        all_before_dict[rep_class] = list(all_before_df_win[0])
        after_dict[rep_class] = list(after_df_win[0])
        all_after_dict[rep_class] = list(all_after_df_win[0])
        
        # make polynomial features
        polynomial_features= PolynomialFeatures(degree = 4)
        ## I tried a few degrees but 4 looks the most reasonable
        x_before_reshape = np.array(x_before).reshape(-1,1)
        # print(f"after reshaping: {x_before_reshape.shape}")
        p_before = polynomial_features.fit_transform(x_before_reshape)
        x_after_reshape = np.array(x_after).reshape(-1,1)
        # print(f"after reshaping: {x_after_reshape.shape}")
        p_after = polynomial_features.fit_transform(x_after_reshape)

        # model polynomial regression
        before_model = sm.WLS(before_dict[rep_class], p_before, weights=np.full_like(before_dict[rep_class], num_sig_transcripts)).fit()
        all_before_model = sm.WLS(all_before_dict[rep_class], p_before, weights=np.full_like(all_before_dict[rep_class], num_all_transcripts)).fit()
        after_model = sm.WLS(after_dict[rep_class], p_after, weights=np.full_like(after_dict[rep_class], num_sig_transcripts)).fit()
        all_after_model = sm.WLS(all_after_dict[rep_class], p_after, weights=np.full_like(all_after_dict[rep_class], num_all_transcripts)).fit()

        ## calculate predicted polynomial
        before_ypred = before_model.predict(p_before) 
        all_before_ypred = all_before_model.predict(p_before) 
        after_ypred = after_model.predict(p_after) 
        all_after_ypred = all_after_model.predict(p_after) 

        ## calculate predicted confidence interval
        _, upper_before_model,lower_before_model = wls_prediction_std(before_model, alpha=0.05)
        _, upper_all_before_model,lower_all_before_model = wls_prediction_std(all_before_model, alpha=0.05)
        _, upper_after_model,lower_after_model = wls_prediction_std(after_model, alpha=0.05)
        _, upper_all_after_model,lower_all_after_model = wls_prediction_std(all_after_model, alpha=0.05)

        rep_label = rep_class.replace("_", " ")
        ax.plot(x_before_reshape, before_dict[rep_class], color = colors[rep_class], linewidth=2)
        ax.plot(x_after_reshape, after_dict[rep_class], color = colors[rep_class], linewidth=2)
        ax.plot(x_before_reshape, before_ypred, label = f"{rep_label} pol. reg.\nforeground", color = colors[rep_class], linewidth = 3, linestyle = (0, (5, 2)))
        ax.plot(x_after_reshape, after_ypred, color = colors[rep_class], linewidth=3, linestyle = (0, (5, 2)))
        ax.fill_between(x_before, upper_before_model,lower_before_model, color=colors[rep_class], alpha = 0.3)
        ax.fill_between(x_after, upper_after_model,lower_after_model, color=colors[rep_class], alpha = 0.3)
        
        if all_before_dict !={} and all_after_dict !={} and all_transcripts!=0:
            max_before = max(all_before_dict[rep_class])
            max_after = max(all_after_dict[rep_class])
            if max_before>max_percentage:
                max_percentage=max_before
            if max_after>max_percentage:
                max_percentage=max_after

            ax.plot(x_before, all_before_dict[rep_class], color = colors[rep_class], linestyle = (0, (1, 2)), linewidth = 2)                    
            ax.plot(x_after, all_after_dict[rep_class], color = colors[rep_class], linestyle = (0, (1, 2)), linewidth = 2)   
            ax.plot(x_before, all_before_ypred, color = colors[rep_class],label = f"{rep_label} pol. reg.\nbackground", linewidth=1, linestyle = (0, (5, 2)))
            ax.plot(x_after, all_after_ypred, color = colors[rep_class], linewidth=1, linestyle = (0, (5, 2)))
            ax.fill_between(x_before, upper_all_before_model,lower_all_before_model, color=colors[rep_class], alpha = 0.2)
            ax.fill_between(x_after, upper_all_after_model,lower_all_after_model, color=colors[rep_class], alpha = 0.2, label = "confidence\ninterval")           

        if True: ## print model output
            o = sys.stdout
            with open(f'{modelstats_filename}_{rep_class}.txt', 'w') as f:
                sys.stdout = f
                print(f"Species: {species} \nsummary of statistical models from polynomial regression for repeat class: {rep_label}")
                print(f"\n############## before transcript foreground model")
                print(before_model.summary())
                print(f"\n\n############## before transcript background model")
                print(all_before_model.summary())
                print(f"\n\n############## after transcript foreground model")
                print(after_model.summary())
                print(f"\n\n############## after transcript background model")
                print(all_after_model.summary())
            sys.stdout = o

        max_percentage=int(max_percentage*1.3)
        if max_percentage == 0 or max_percentage>100:
            max_percentage= 100

        plt.vlines(x= 0, ymin=0, ymax=max_percentage, colors="#000000", linestyles="dashed", label="transcript border", linewidth=3)
        plt.xticks(range(-num_bp, num_bp+1, int(num_bp/5)), fontsize = fs)
        plt.yticks(range(0, max_percentage+1, 10), fontsize = fs)

        if legend:
            # plot color legend
            ax.set_xlim([-num_bp, num_bp*1.55])
            legend_colors = ax.legend(loc = "center right", fontsize = fs)
            plt.gca().add_artist(legend_colors)
            # plt.title(f"{species} transcript surroundings {num_bp} bp up and downstream \n({num_sig_transcripts} significant transcripts of {all_transcripts} in CAFE analysis)", fontsize = fs*1.25)
            
        # plot dotted/bold legend
        solid = Line2D([0], [0], color='black', linestyle='-', linewidth=2)
        dotted = Line2D([0], [0], color='black', linestyle=':', linewidth=2)
        handles = [solid, dotted]
        labels = []
        
        if True:
            labels.append(f"foreground transcripts ({num_sig_transcripts})")
            labels.append(f"background transcripts ({all_transcripts})")
        else:
            labels.append(f"significant transcripts ({num_sig_transcripts})")
            labels.append(f"all CAFE transcripts ({all_transcripts})")
        plt.legend(handles, labels, loc = "upper left", fontsize = fs, title=legend_title, title_fontsize=fs)

        plt.title(f"{species} transcript surroundings {num_bp} bp up and downstream\nrepeat category: {rep_label} with polynomial regression and 95% confidence interval", fontsize = fs*1.25)
        plt.xlabel(f"basepairs upstream and downstream from transcript", fontsize = fs)

        plt.ylabel(f"percent of transcripts in which this base is a repeat", fontsize = fs)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 99 or x<1 else f'{int(x)}%'))
        
        plt.tight_layout()
        filename_class = f"{filename}_{rep_class}.png"
        plt.savefig(filename_class, dpi = 300, transparent = plot_transparent_bg)
        plt.close(fig)
        print(f"\t * repeat category {rep_label}    \t--> Figure saved as: {filename_class}")


