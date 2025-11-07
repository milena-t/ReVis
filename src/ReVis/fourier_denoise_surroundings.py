"""
calculate and plot the confidence intervals of the individual TE categories in the transcript surroundings
"""


import parse_gff as gff
from statistics_surroundings import manual_yticks

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D

from sklearn.preprocessing import PolynomialFeatures
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std


def make_fft_of_rep_prop(rep_class_dict, psd_max = 100):
    n = len(rep_class_dict)
    fhat = np.fft.fft(rep_class_dict, n)
    PSD_fft = fhat * np.conj(fhat) / n # Power spectral density
    freq = (1/(n))*np.arange(n) # vector of all frequencies
    L_ind = np.arange(1, np.floor(n/2), dtype = 'int') # indices

    indices = PSD_fft > psd_max # indices 
    PSDClean = PSD_fft*indices 
    fhat = indices * fhat
    ffilt = np.fft.ifft(fhat)

    ### test imaginary component --> should only be noise with like 1e-13 or lower
    if False:
        ffilt_test = ffilt
        print(np.max(np.abs(ffilt_test.imag)))

    ffilt = ffilt.real
    return [PSD_fft,freq,L_ind], ffilt

def plot_confidence_intervals(before_filepath:str, after_filepath:str, num_sig_transcripts:int, num_all_transcripts:int, general_legend_names = True, all_before_filepath:str = "", all_after_filepath:str = "", filename = "cumulative_repeat_presence_around_transcripts_95_perc_CI", modelstats_filename = "pol_reg_sum.txt", plot_fourier_transform = True, legend = True, plot_white_bg = False):
    """
    make a separate plot for each TE category with the foreground/background and before/after lines in it 
    do the polynomial regressions and see the 95% confidence intervals
    I followed this tutorial: https://ostwalprasad.github.io/machine-learning/Polynomial-Regression-using-statsmodel.html
    """
    plot_transparent_bg = True
    if plot_white_bg:
        plot_transparent_bg = False

    out_dir_name = "/".join(filename.split("/")[:-1])
    print(f"individual repeat category plots saved in dir: {out_dir_name}")
    
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
        all_before_dict = { key : [int(v)/num_all_transcripts*100 for v in value] for key, value in all_before_dict.items()}
        all_after_dict = gff.read_dict_from_file(all_after_filepath)
        all_after_dict = { key : [int(v)/num_all_transcripts*100 for v in value] for key, value in all_after_dict.items()}
        
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
    x_before = range(-num_bp, 0)
    x_after = range(num_bp)

    for rep_class in rep_classes:

        rep_label = rep_class.replace("_", " ")
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

        ####### POLYNOMIAL REGRESSION
        
        ### fourier-denoise the data
        psd_max = 50
        before_model_psd_plot_list, before_fourier_denoise = make_fft_of_rep_prop(before_dict[rep_class], psd_max = psd_max)
        all_before_model_psd_plot_list, all_before_fourier_denoise = make_fft_of_rep_prop(all_before_dict[rep_class], psd_max = psd_max)
        after_model_psd_plot_list, after_fourier_denoise = make_fft_of_rep_prop(after_dict[rep_class], psd_max = psd_max)
        all_after_model_psd_plot_list, all_after_fourier_denoise = make_fft_of_rep_prop(all_after_dict[rep_class], psd_max = psd_max)


        # make polynomial features
        ## I tried a few degrees but 4 looks the most reasonable
        polynomial_features= PolynomialFeatures(degree = 4)

        x_before_reshape = np.array(x_before).reshape(-1,1)
        p_before = polynomial_features.fit_transform(x_before_reshape)
        x_after_reshape = np.array(x_after).reshape(-1,1)
        p_after = polynomial_features.fit_transform(x_after_reshape)

        # model polynomial regression
        before_model = sm.WLS(list(before_fourier_denoise), p_before, weights=np.full_like(before_dict[rep_class], num_sig_transcripts)).fit()
        all_before_model = sm.WLS(list(all_before_fourier_denoise), p_before, weights=np.full_like(all_before_dict[rep_class], num_all_transcripts)).fit()
        after_model = sm.WLS(list(after_fourier_denoise), p_after, weights=np.full_like(after_dict[rep_class], num_sig_transcripts)).fit()
        all_after_model = sm.WLS(list(all_after_fourier_denoise), p_after, weights=np.full_like(all_after_dict[rep_class], num_all_transcripts)).fit()

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

        ax.plot(x_before_reshape, before_fourier_denoise, color = colors[rep_class], linewidth=2)
        ax.plot(x_after_reshape, after_fourier_denoise, color = colors[rep_class], linewidth=2)
        ax.plot(x_before_reshape, before_ypred, label = f"{rep_label} pol. reg.\nforeground", color = colors[rep_class], linewidth = 3, linestyle = (0, (5, 2)))
        ax.plot(x_after_reshape, after_ypred, color = colors[rep_class], linewidth=3, linestyle = (0, (5, 2)))
        ax.fill_between(x_before, upper_before_model,lower_before_model, color=colors[rep_class], alpha = 0.3)
        ax.fill_between(x_after, upper_after_model,lower_after_model, color=colors[rep_class], alpha = 0.3)
        
        if all_before_dict !={} and all_after_dict !={} and num_all_transcripts!=0:
            max_before = max(all_before_dict[rep_class])
            max_after = max(all_after_dict[rep_class])
            if max_before>max_percentage:
                max_percentage=max_before
            if max_after>max_percentage:
                max_percentage=max_after

            ax.plot(x_before, all_before_fourier_denoise, color = colors[rep_class], linestyle = (0, (1, 3)), linewidth = 2)                    
            ax.plot(x_after, all_after_fourier_denoise, color = colors[rep_class], linestyle = (0, (1, 3)), linewidth = 2)   
            ax.plot(x_before, all_before_ypred, color = colors[rep_class],label = f"{rep_label} pol. reg.\nbackground", linewidth=1, linestyle = (0, (10, 3)))
            ax.plot(x_after, all_after_ypred, color = colors[rep_class], linewidth=1, linestyle = (0, (10, 3)))
            ax.fill_between(x_before, upper_all_before_model,lower_all_before_model, color=colors[rep_class], alpha = 0.2)
            ax.fill_between(x_after, upper_all_after_model,lower_all_after_model, color=colors[rep_class], alpha = 0.2, label = "confidence\ninterval")           

        if False: ## print model output
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

        filename_out = f"{filename}_{rep_class}.png"
        
        max_percentage=int(max_percentage*1.3)
        if max_percentage == 0 or max_percentage>100:
            max_percentage= 100

        plt.vlines(x= 0, ymin=0, ymax=max_percentage, colors="#000000", linestyles="dashed", label="transcript border", linewidth=3)
        plt.xticks(range(-num_bp, num_bp+1, int(num_bp/5)), fontsize = fs)
        manual_yticks(max_percentage, ax, fs)

        if legend:
            # plot color legend
            ax.set_xlim([-num_bp, num_bp*1.55])
            legend_colors = ax.legend(loc = "center right", fontsize = fs)
            plt.gca().add_artist(legend_colors)
            
        # plot dotted/bold legend
        solid = Line2D([0], [0], color=colors[rep_class], linestyle='-', linewidth=2)
        dotted = Line2D([0], [0], color=colors[rep_class], linestyle=':', linewidth=2)
        handles = [solid, dotted]
        labels = []
        if True:
            labels.append(f"foreground transcripts ({num_sig_transcripts})")
            labels.append(f"background transcripts ({num_all_transcripts})")
        else:
            labels.append(f"significant transcripts ({num_sig_transcripts})")
            labels.append(f"all CAFE transcripts ({num_all_transcripts})")
        plt.legend(handles, labels, loc = "upper left", fontsize = fs, title = "fourier denoised", title_fontsize = fs)

        plt.title(f"{species} transcript surroundings {num_bp} bp up and downstream\nrepeat category: {rep_label} with polynomial regression and 95% confidence interval", fontsize = fs*1.25)
        plt.xlabel(f"basepairs upstream and downstream from transcript", fontsize = fs)

        plt.ylabel(f"percent of transcripts in which this base is a repeat", fontsize = fs)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 99 or x<1 else f'{int(x)}%'))
        
        plt.tight_layout()
        filename_out = f"{filename}_95perc_CI_fft_{rep_class}.png"
        plt.savefig(filename_out, dpi = 300, transparent = plot_transparent_bg)
        plt.close(fig)
        filename_out = filename_out.split("/")[-1]
        print(f"\t * repeat category {rep_label} \t--> Figure saved as: {filename_out}")

        if plot_fourier_transform:

            max_percentage = 0
            if legend:
                fig, ax = plt.subplots(1, 1, figsize=(23, 10))
            else:
                fig, ax = plt.subplots(1, 1, figsize=(21, 12))

            if False:
                ### plot the power spectral density function to find the upper limit for the fourier denoising

                def plot_ft_power_spectral_density(in_list, upper_lim = psd_max):
                    """
                    input is a list of [before_model, after_model, all_before_model, all_after_model] where each element is 
                    again a list of [PSD_fft,freq,L_ind] returned by make_fft_of_rep_prop
                    """
                    fig, ax = plt.subplots(1,1, figsize=(15, 5))
                    colors = ["#c4bbaf", "#a5978b", "#5c4742", "#8d5b4c",]
                    labels = ["before_model", "after_model", "all_before_model", "all_after_model"]

                    for i, section in enumerate(in_list):    
                        PSD_fft,freq,L_ind = section
                        ax.plot(freq[L_ind], PSD_fft[L_ind], color = colors[i], linewidth = 2, label = labels[i])
                        ax.set_xlim([freq[L_ind[0]], freq[L_ind[-1]]])
                        ax.set_ylim([0, upper_lim])
                    plt.title(f"{rep_class}")
                    plt.legend()
                    filename_PSD = filename.replace("_95_perc_CI", "_psd")
                    filename_PSD = f"{filename_PSD}_{rep_class}.png"
                    plt.savefig(filename_PSD, dpi = 300, transparent = plot_transparent_bg)
                    print(f"{rep_class} PSD plot saved as: {filename_PSD}")

                plot_ft_power_spectral_density([before_model_psd_plot_list, after_model_psd_plot_list, all_before_model_psd_plot_list, all_after_model_psd_plot_list])
                # raise RuntimeError(f"test until {rep_class}")
            
            raw_transparency = 0.4
            ax.plot(x_before, before_fourier_denoise, label = f"fourier\ndenoise\nforeground", color = colors[rep_class], linewidth = 3)#, linestyle = (0, (5, 2)))
            ax.plot(x_before, before_dict[rep_class], color = colors[rep_class], linewidth=2, alpha = raw_transparency)
            ax.plot(x_after, after_fourier_denoise, color = colors[rep_class], linewidth=3)#, linestyle = (0, (5, 2)))
            ax.plot(x_after, after_dict[rep_class], color = colors[rep_class], linewidth=2, alpha = raw_transparency)

            if all_before_dict !={} and all_after_dict !={} and num_all_transcripts!=0:
                max_before = max(all_before_dict[rep_class])
                max_after = max(all_after_dict[rep_class])
                if max_before>max_percentage:
                    max_percentage=max_before
                if max_after>max_percentage:
                    max_percentage=max_after

                ax.plot(x_before, all_before_fourier_denoise, color = colors[rep_class],label = f"fourier\ndenoise\nbackground", linewidth=1 )#, linestyle = (0, (5, 2)))
                ax.plot(x_before, all_before_dict[rep_class], color = colors[rep_class], linestyle = (0, (1, 10)), linewidth = 2, alpha = raw_transparency)                    
                ax.plot(x_after, all_after_fourier_denoise, color = colors[rep_class], linewidth=1 )#, linestyle = (0, (5, 2)))
                ax.plot(x_after, all_after_dict[rep_class], color = colors[rep_class], linestyle = (0, (1, 10)), linewidth = 2, alpha = raw_transparency)   


            filename_fourier_denoise = f"{filename}_fft_denoise_{rep_class}.png"

            max_percentage=int(max_percentage*1.3)
            if max_percentage == 0 or max_percentage>100:
                max_percentage= 100

            plt.vlines(x= 0, ymin=0, ymax=max_percentage*2, colors="#000000", linestyles="dashed", label="transcript border", linewidth=3)
            plt.xticks(range(-num_bp, num_bp+1, int(num_bp/5)), fontsize = fs)
            
            manual_yticks(max_percentage, ax, fs)
            
            ax.set_ylim([0,max_percentage*1.5])

            if legend:
                # plot color legend
                ax.set_xlim([-num_bp, num_bp*1.55])
                legend_colors = ax.legend(loc = "center right", fontsize = fs)
                plt.gca().add_artist(legend_colors)
                
            # plot dotted/bold legend
            solid = Line2D([0], [0], color=colors[rep_class], linestyle='-', linewidth=2, alpha = raw_transparency)
            dotted = Line2D([0], [0], color=colors[rep_class], linestyle=':', linewidth=2, alpha = raw_transparency)
            handles = [solid, dotted]
            labels = []
            if True:
                labels.append(f"foreground transcripts ({num_sig_transcripts})")
                labels.append(f"background transcripts ({num_all_transcripts})")
            else:
                labels.append(f"significant transcripts ({num_sig_transcripts})")
                labels.append(f"all CAFE transcripts ({num_all_transcripts})")
            plt.legend(handles, labels, loc = "upper left", fontsize = fs)

            plt.title(f"{species} transcript surroundings {num_bp} bp up and downstream\nrepeat category: {rep_label} with polynomial regression and 95% confidence interval", fontsize = fs*1.25)
            plt.xlabel(f"basepairs upstream and downstream from transcript", fontsize = fs)

            plt.ylabel(f"percent of transcripts in which this base is a repeat", fontsize = fs)
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 99 or x<1 else f'{int(x)}%'))
            
            plt.tight_layout()
            # filename_out = f"{filename}_{rep_class}.png"
            plt.savefig(filename_fourier_denoise, dpi = 300, transparent = plot_transparent_bg)
            plt.close(fig)
            filename_fourier_denoise = filename_fourier_denoise.split("/")[-1]
            print(f"\t\t\t\t\t--> Fourier figure saved: {filename_fourier_denoise}")


