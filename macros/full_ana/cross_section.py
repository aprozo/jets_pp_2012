#!/usr/bin/env python
# coding: utf-8

# In[12]:


import ROOT
get_ipython().run_line_magic('jsroot', 'on')
# calculate events ratio from vertex histograms
path='/home/prozorov/dev/star/jets_pp_2012/output/'
workdir='/home/prozorov/dev/star/jets_pp_2012/macros/full_ana/'

triggers = ['JP2', 'HT2']
jetR = [0.2, 0.3, 0.4, 0.5, 0.6]
jetR='0.5'


bad_runs_JP2=[  13048019, 13048092, 13048093,
                13049006, 13049007, 13051074,
                13051099, 13052061, 13061035,
                13064067, 13068060, 13069023,
                13070061, 13063012]


# In[13]:


# print ROOT version
print(ROOT.gROOT.GetVersion())


# ```bash
# export PATH=/home/prozorov/install/RooUnfold:${PATH}
# export PYTHONPATH=/home/prozorov/install/RooUnfold:${PYTHONPATH}
# export LD_LIBRARY_PATH=/home/prozorov/install/RooUnfold:${LD_LIBRARY_PATH}
# ```
# 

# In[14]:


ROOT.gSystem.Load("/home/prozorov/install/RooUnfold/lib/libRooUnfold")  # load RooUnfold library


# In[15]:


from ROOT import RooUnfoldResponse


# ### plot Vertex Z distributions

# In[16]:


import math
# calculate events ratio from vertex histograms
event_ratio = {}
ROOT.gStyle.SetOptStat(100)
canvas = ROOT.TCanvas('canvas', '', 1200, 600)
canvas.Divide(3, 1)
canvas.Draw()

line = ROOT.TLine()
line.SetLineColor(ROOT.kRed)
line.SetLineStyle(2)

text =ROOT.TLatex()
text.SetTextSize(0.03)
text.SetTextColor(ROOT.kRed)

for i, trigger in enumerate(triggers):
    file = ROOT.TFile(f'/home/prozorov/dev/star/jets_pp_2012/output/merged_data_{trigger}_R{jetR}.root')
    vertex_hist = file.Get('QA_histograms/vz')
    vertex_hist.RebinX()
    vertex_hist.SetName(f'hVertexZ_{trigger}')
    vertex_hist.GetXaxis().SetTitle('Z_{vertex} [cm]')
    vertex_hist.SetDirectory(0)
    vertex_hist.SetTitle(trigger)
    vertex_hist.GetXaxis().SetRangeUser(-60, 60)
    vertex_hist.GetYaxis().SetRangeUser(0 ,vertex_hist.GetMaximum() * 1.05)

    canvas.cd(i + 1)
    pad=ROOT.gPad
    ROOT.gPad.SetLeftMargin(0.05)
    ROOT.gPad.SetRightMargin(0.01)
    vertex_hist.Draw('hist')
    gaus_fit = ROOT.TF1(f"gaus_{trigger}", "gaus", -30, 30)
    gaus_fit.SetLineColor(ROOT.kPink)
    vertex_hist.Fit(gaus_fit, "REMQ")
    gaus_fit.Draw('same')

    binwidth = vertex_hist.GetBinWidth(2)
    integral = vertex_hist.Integral(vertex_hist.FindBin(-30), vertex_hist.FindBin(30)) * binwidth
    total_integral = math.sqrt(2 * math.pi) * gaus_fit.GetParameter(0) * gaus_fit.GetParameter(2)
    ratio = integral / total_integral
    event_ratio[trigger] = ratio

    line.DrawLine(-30, 0, -30, 1.0 * vertex_hist.GetMaximum())
    line.DrawLine(30, 0, 30, 1.0 * vertex_hist.GetMaximum())
    text.DrawLatexNDC(0.45, 0.7, f'{100*ratio:.1f} % total')


canvas.SaveAs('plots/vertex_ratio.pdf')


# In[19]:


import ROOT
import uproot
import awkward as ak
import numpy as np
path='/home/prozorov/dev/star/jets_pp_2012/output/'
triggers = ['JP2', 'HT2']
tree_name = 'ResultTree'  # Update this to your tree name
branch_names = ['runid1', 'pt', 'trigger_match_JP2', 'trigger_match_HT2']
pt_vs_run_raw = {}
pt_vs_run_trigger = {}
# read run_map from file
run_map = {}
run_map_inverse = {}

with open('/home/prozorov/dev/star/jets_pp_2012/macros/ana/run_map.txt', 'r') as f:
    for line in f:
        run, i = line.strip().split()
        run_map[int(run)] = int(i)
        run_map_inverse[int(i)] = int(run)

import json
with open("my_parameters.json") as f:
    pt_reco_bins = json.load(f)["pt_reco_bins"]

pt_reco_bins=np.array(pt_reco_bins)

for trigger in triggers:
    print(f'Processing trigger {trigger}')
    file_name = path + f'merged_data_{trigger}_R{jetR}.root'
    my_file=ROOT.TFile(file_name)
    analyzed_event_hist= my_file.Get('hEventsRun')
    nbins=analyzed_event_hist.GetNbinsX()
    root_file=uproot.open(file_name)
    tree = root_file[tree_name]
    tree_arrays = tree.arrays(branch_names)

    # convert runid1 to run_mapped
    tree_arrays['run_mapped'] = [run_map[runid] for runid in tree_arrays['runid1']]

    if trigger in ['JP2', 'HT2']:
        df = ak.to_rdataframe({"pt": tree_arrays['pt'], "run_mapped": tree_arrays['run_mapped'], "trigger_match": tree_arrays[f'trigger_match_{trigger}']})
    else:
        df = ak.to_rdataframe({"pt": tree_arrays['pt'], "run_mapped": tree_arrays['run_mapped']})


    title = '; run ;p_{t}, GeV/c; counts'

    pt_vs_run = df.Histo2D((f'{trigger}_pt_vs_run', f'{trigger} pt'+title,nbins,0,nbins, len(pt_reco_bins)-1, pt_reco_bins), "run_mapped","pt")


        # set labels equal to runid from run_map
    # for i in range(1, pt_vs_run.GetNbinsX()+1):
    #     pt_vs_run.GetXaxis().SetBinLabel(i, str(run_map_inverse[i-1]))

    pt_vs_run_raw[trigger] = pt_vs_run

    if trigger in ['JP2', 'HT2']:
        hist_2d_trigger= df.Define("triggered_pt", "pt[trigger_match]").Histo2D((f'{trigger}_pt_vs_run_trigger', f'{trigger} pt triggered'+title,nbins,0,nbins, len(pt_reco_bins)-1, pt_reco_bins), "run_mapped","triggered_pt")

        pt_vs_run_trigger[trigger] = hist_2d_trigger


outfile_2d = ROOT.TFile(workdir+"pt_vs_run_binning.root", "RECREATE")
for trigger in triggers:
    pt_vs_run_raw[trigger].Write()
    if trigger in ['JP2', 'HT2']:
        pt_vs_run_trigger[trigger].Write()
outfile_2d.Close()


# ### Read raw pt distributions

# In[21]:


pt_vs_run_raw = {}
pt_vs_run_trigger = {}

infile_2d = ROOT.TFile(workdir+"pt_vs_run_binning.root", "READ")
for trigger in triggers:
    pt_vs_run_raw[trigger] = infile_2d.Get(f'{trigger}_pt_vs_run')
    pt_vs_run_raw[trigger].SetDirectory(0)
    if trigger in ['JP2', 'HT2']:
        pt_vs_run_trigger[trigger] = infile_2d.Get(f'{trigger}_pt_vs_run_trigger')
        pt_vs_run_trigger[trigger].SetDirectory(0)

# pt_vs_run_trigger['MB']= pt_vs_run_raw['MB']


# ### Get analyzed events histograms

# In[24]:


analyzed_events = {}
for trigger in triggers:
    file_name = path + f'merged_data_{trigger}_R{jetR}.root'

    my_file=ROOT.TFile(file_name)
    analyzed_events[trigger]= my_file.Get('hEventsRun')
    analyzed_events[trigger].SetDirectory(0)


# ### Read all prescales, livetimes, maxnevents, luminosities

# In[ ]:


# get 1d histograms - presacles, livetime, maxnevents
lumi=ROOT.TFile(workdir+"lumi.root")
prescale={}
livetime={}
maxnevent={}
luminosity={}
for trigger in triggers:
    prescale[trigger]= lumi.Get(f'prescale_{trigger}')
    livetime[trigger]= lumi.Get(f'livetime_{trigger}')
    maxnevent[trigger]= lumi.Get(f'nevents_{trigger}')
    luminosity[trigger]= lumi.Get(f'luminosity_{trigger}')


# ### Get Trigger efficiency from embedding

# In[ ]:


trigger_efficiency={}
trigger_mc_efficiency={}

# load trigger efficiency from embedding
for trigger in ['JP2', 'HT2']:
    input_file = ROOT.TFile(workdir+f"/plots/embedding_root/{trigger}_over_reconstructed_reco.root", "READ")
    trigger_efficiency[trigger]=input_file.Get(f"h_ratio").Clone(f'h_ratio_{trigger}_reco')
    trigger_efficiency[trigger].SetDirectory(0)

    input_file_mc = ROOT.TFile(workdir+f"/plots/embedding_root/{trigger}_over_reconstructed_mc.root", "READ")
    trigger_mc_efficiency[trigger]=input_file_mc.Get(f"h_ratio").Clone(f'h_ratio_{trigger}_mc')
    trigger_mc_efficiency[trigger].SetDirectory(0)

input_file=ROOT.TFile(workdir+f"/plots/embedding_root/reconstructed_over_all_reco.root", "READ")
reconstruction_efficiency=input_file.Get(f"h_ratio").Clone(f'reconstruction_efficiency')
reconstruction_efficiency.SetDirectory(0)

input_file_mc=ROOT.TFile(workdir+f"/plots/embedding_root/reconstructed_over_all_mc.root", "READ")
reconstruction_mc_efficiency=input_file_mc.Get(f"h_ratio").Clone(f'reconstruction_mc_efficiency')
reconstruction_mc_efficiency.SetDirectory(0)


print("Trigger efficiency loaded")


# In[ ]:


import math
def divide(num, denom):
    hnew = num.Clone()  # Clone the first histogram
    for j in range(1, hnew.GetNbinsY() + 1):
        for i in range(1, hnew.GetNbinsX() + 1):
            bin_content = denom.GetBinContent(i)
            if bin_content == 0:
                scale = 0
            else:
                scale = 1 / bin_content
            hnew.SetBinContent(i, j, hnew.GetBinContent(i, j) * scale)
            hnew.SetBinError(i, j, hnew.GetBinError(i, j) * scale)
    hnew.SetDirectory(0)  # Detach from any file
    return hnew

def multiply(num, denom):
    hnew = num.Clone()  # Clone the first histogram
    for j in range(1, hnew.GetNbinsY() + 1):
        for i in range(1, hnew.GetNbinsX() + 1):
            scale = denom.GetBinContent(i)
            hnew.SetBinContent(i, j, hnew.GetBinContent(i, j) * scale)
            hnew.SetBinError(i, j, hnew.GetBinError(i, j) * scale)
    hnew.SetDirectory(0)  # Detach from any file
    return hnew

# make 1d histogram out of this histograms
def make_run_projection(hist):
    counter_runs=0
    name = hist.GetName()
    temp = hist.ProjectionY('', 1, 1)
    output=temp.Clone(f'{name}')
    for i in range(1, hist.GetNbinsX() + 1):
         projection = hist.ProjectionY(f'{name}_{i}', i, i)
         runid = int(hist.GetXaxis().GetBinLabel(i))
         # check if run is in bad runs
         if runid in bad_runs_JP2:
             continue
         if projection.Integral() == 0:
            continue
         counter_runs+=1
         output.Add(projection)
    output.Scale(1/counter_runs)
    output.SetDirectory(0)  # Detach from any file

    return output

# go by bin by bin and multiply by trigger efficiency interpolated
def divide_by_efficiency(data, eff):
    for i in range(1, data.GetNbinsX()+1):
        bin_content = data.GetBinContent(i)
        bin_error = data.GetBinError(i)
        bin_center = data.GetBinCenter(i)
        # get trigger efficiency
        eff_value = eff.GetBinContent(eff.GetXaxis().FindBin(bin_center))
        eff_error = eff.GetBinError(eff.GetXaxis().FindBin(bin_center))
        # check if trigger efficiency is 0
        if eff_value == 0:
            eff_value = 1
        # set new value
        data.SetBinContent(i, bin_content /eff_value)
        # calculate uncertainty propagation
        # check if not 0
        if bin_content == 0:
            continue
        new_bin_error = bin_content / eff_value * math.sqrt((bin_error/bin_content)**2 + (eff_error/eff_value)**2)
        data.SetBinError(i, new_bin_error)
    return data

def multiply_by_efficiency (data,eff):
    for i in range(1, data.GetNbinsX()+1):
        bin_content = data.GetBinContent(i)
        bin_error = data.GetBinError(i)
        bin_center = data.GetBinCenter(i)
        # get trigger efficiency
        eff_value = eff.GetBinContent(eff.GetXaxis().FindBin(bin_center))
        eff_error = eff.GetBinError(eff.GetXaxis().FindBin(bin_center))
        # check if trigger efficiency is 0
        if eff_value == 0:
            eff_value = 1
        # set new value
        data.SetBinContent(i, bin_content *eff_value)
        # calculate uncertainty propagation
        # check if not 0
        if bin_content == 0:
            continue
        new_bin_error = bin_content * eff_value * math.sqrt((bin_error/bin_content)**2 + (eff_error/eff_value)**2)
        data.SetBinError(i, new_bin_error)
    return data

def divide_by_binwidth(data):
    for i in range(1, data.GetNbinsX()+1):
        bin_content = data.GetBinContent(i)
        bin_error = data.GetBinError(i)
        bin_width = data.GetBinWidth(i)
        # set new value
        data.SetBinContent(i, bin_content /bin_width)
        data.SetBinError(i, bin_error /bin_width)
    return data


# ## Unfolding

# In[ ]:


# responseFile=ROOT.TFile(f"/home/prozorov/dev/star/jets_pp_2012/macros/unfolding/response_{trigger}.root")
responses={}

for trigger in triggers:
        responseFile=ROOT.TFile(f"/home/prozorov/dev/star/jets_pp_2012/macros/unfolding/root/response{trigger}R{jetR}.root")
        responses[trigger] = responseFile.Get("my_response")
        responses[trigger].SetName(f"response_{trigger}")


rooUnfoldResponse = responses['MB']                                      


# In[ ]:


from ROOT import RooUnfoldBayes
temp1d = pt_vs_run_trigger['JP2'].ProjectionY()
# temp1d=responseFile.Get("MeasuredTest")
temp1d=divide_by_efficiency (temp1d, trigger_efficiency['JP2'])
temp1d.SetLineColor(ROOT.kBlue)
unfolding=ROOT.RooUnfoldBayes(rooUnfoldResponse, temp1d, 5)
unfolded = unfolding.Hunfold()
unfolded.SetLineColor(ROOT.kRed)
# divide by bin width
unfolded=divide_by_binwidth(unfolded)
temp1d=divide_by_binwidth(temp1d)


# In[ ]:


can=ROOT.TCanvas('can', '', 1000, 600)
can.Draw()
can.SetLogy()
temp1d.SetTitle('JP2 trigger')
temp1d.GetYaxis().SetTitle('dN/dp_{t}')
temp1d.GetYaxis().SetTitleOffset(0.8)
temp1d.Draw('hist')
unfolded.Draw('hist same')
leg=ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(temp1d, 'raw data', 'l')
leg.AddEntry(unfolded, 'unfolded', 'l')
leg.Draw()
can.SaveAs('plots/unfolded.pdf')


# In[ ]:


# print x-axis edges vector 
xaxis = temp1d.GetXaxis()
x_edges = []
for i in range(1, temp1d.GetNbinsX() + 2):
    x_edges.append(xaxis.GetBinLowEdge(i))
print(x_edges)


# In[ ]:


# load Dmitry's histograms
file=ROOT.TFile(workdir+'./jet_cross_section_dmitriy.root')
base_hist=file.Get('crossSection_systematic')
# base_hist=file.Get('crossSection_statistic')
base_hist.SetDirectory(0)
def divide_different_bins (num, denom):
    hnew = num.Clone()  # Clone the first histogram
    for i in range(1, hnew.GetNbinsX() + 1):
        scale = denom.Interpolate(hnew.GetXaxis().GetBinCenter(i))
        # scale =denom.GetBinContent( denom.GetXaxis().FindBin(hnew.GetXaxis().GetBinCenter(i)))
        if scale == 0:
            scale = 1
        hnew.SetBinContent(i, hnew.GetBinContent(i) /scale)
        hnew.SetBinError(i,  hnew.GetBinError(i) / scale)
    hnew.SetDirectory(0)  # Detach from any file
    return hnew


# ### Combine all histograms on per-run basis

# In[ ]:


my_cross_section = {}
my_cross_section_no_unfolding = {}

for trigger in triggers:
    print(f'Calculating cross section for {trigger}')

    temp = pt_vs_run_trigger[trigger]
    # =============================================
    ## normalize by analyzed events
    temp = divide(temp, analyzed_events[trigger]) 

    #  by maxnevent
    temp = multiply(temp, maxnevent[trigger]) 

    # by luminosity
    temp = divide(temp, luminosity[trigger]) 

    # scale by ratio of passed events
    temp.Scale(1./event_ratio[trigger])

    # make run projection
    temp = make_run_projection(temp)
   # scale by acceptance pseudo-rapidity
    temp.Scale(1./(2*0.4))

    # trigger efficiency
    if trigger in ['JP2', 'HT2']:
        temp = divide_by_efficiency(temp, trigger_efficiency[trigger])

    # temp = divide_by_efficiency(temp, reconstruction_efficiency)


    # make unfolding
    temp_nounfolding=temp.Clone(f'{trigger}_cross_section_no_unfolding')
    temp=temp_nounfolding.Clone(f'{trigger}_cross_section')
    # # # unfolding = ROOT.RooUnfoldBayes(responses[trigger], temp, 5 )
    # unfolding = ROOT.RooUnfoldBayes(rooUnfoldResponse, temp, 5 )
    # temp = unfolding.Hunfold()
    temp = divide_by_efficiency(temp, reconstruction_mc_efficiency)

    # normalize by pt bin width
    temp = divide_by_binwidth(temp)
    my_cross_section[trigger] = temp

    temp_nounfolding = divide_by_binwidth(temp_nounfolding)
    temp_nounfolding = multiply_by_efficiency(temp_nounfolding, reconstruction_mc_efficiency)
    my_cross_section_no_unfolding[trigger] = temp_nounfolding



# In[ ]:


ROOT.gStyle.SetOptStat(0)
canvas = ROOT.TCanvas('canvasCS', 'canvasCS', 800, 800)
canvas.Divide(1, 2)
# draw all cross sections
colors={
    'JP2': ROOT.kAzure-1,
    'HT2': ROOT.kRed+1,
    'MB': ROOT.kOrange-1,
    'Dmitriy': ROOT.kViolet,
}

# triggers = ['JP2']
triggers = ['HT2', 'JP2', 'MB']

# Upper pad for histograms
pad1 = canvas.cd(1)
pad1.SetPad(0.0, 0.5, 1.0, 1.0)
pad1.SetTopMargin(0.1)
pad1.SetBottomMargin(0.0)
pad1.SetLogy()

canvas.Draw()

canvas.cd(1)
base_hist.GetYaxis().SetRangeUser(1e0+0.2, 1e7)
base_hist.SetMarkerColor(colors['Dmitriy'])
base_hist.SetLineColor(colors['Dmitriy'])
base_hist.SetMarkerStyle(29)
base_hist.SetMarkerSize(2)
base_hist.Draw('E1')

for i, trigger in enumerate(triggers): 
    my_cross_section[trigger].SetLineColor(colors[trigger])
    my_cross_section[trigger].SetMarkerColor( colors[trigger])
    my_cross_section[trigger].SetMarkerStyle(20+i)
    my_cross_section[trigger].Draw('E1 same')

    my_cross_section_no_unfolding[trigger].SetLineColor(colors[trigger]+2)
    my_cross_section_no_unfolding[trigger].SetMarkerColor( colors[trigger]+2)
    my_cross_section_no_unfolding[trigger].SetMarkerStyle(20+i)
    # my_cross_section_no_unfolding[trigger].Draw('E0 same')

legend = ROOT.TLegend(0.7, 0.5, 0.88, 0.88)
for trigger in triggers:
    legend.AddEntry(my_cross_section[trigger], trigger, "lep")
    # legend.AddEntry(my_cross_section_no_unfolding[trigger], f'{trigger} no unfolding', "lep")
legend.AddEntry(base_hist, 'Dmitriy', "lep")
legend.Draw()
#========================
# Lower pad for ratio
pad2 = canvas.cd(2)
pad2.SetPad(0.0, 0.0, 1.0, 0.5)
pad2.SetTopMargin(0.0)
pad2.SetBottomMargin(0.3)

canvas.cd(2)

ratio_hist={}
ratio_hist_no_unfolding={}

for trigger in triggers:
    ratio_hist[trigger] = my_cross_section[trigger].Clone(f"ratio_{trigger}")

    ratio_hist[trigger]=divide_different_bins(ratio_hist[trigger], base_hist)
    ratio_hist[trigger].GetYaxis().SetRangeUser(0, 3.95)
    ratio_hist[trigger].GetXaxis().SetRangeUser(base_hist.GetXaxis().GetXmin()+0.1, base_hist.GetXaxis().GetXmax())
    ratio_hist[trigger].GetYaxis().SetTitle('Ratio to Dmitry')
    ratio_hist[trigger].GetXaxis().SetTitle('p_{t} [GeV/c]')
    ratio_hist[trigger].SetDirectory(0) 
    ratio_hist[trigger].Draw('E1 same')
    ratio_hist_no_unfolding[trigger] = my_cross_section_no_unfolding[trigger].Clone(f"ratio_no_unfolding_{trigger}")
    ratio_hist_no_unfolding[trigger]=divide_different_bins(ratio_hist_no_unfolding[trigger], base_hist)
    ratio_hist_no_unfolding[trigger].GetYaxis().SetRangeUser(0, 3.95)
    ratio_hist_no_unfolding[trigger].GetXaxis().SetRangeUser(base_hist.GetXaxis().GetXmin()+0.1, base_hist.GetXaxis().GetXmax())
    # ratio_hist_no_unfolding[trigger].Draw('E0 same')
     # Detach from any file
# draw line at 1
line = ROOT.TLine()
line.SetLineColor(colors['Dmitriy'])
line.DrawLine(base_hist.GetXaxis().GetXmin(), 1, base_hist.GetXaxis().GetXmax(), 1)
line.Draw()


# In[ ]:


canvas.SaveAs('plots/cross_section_nounfolding.pdf')


# 
# ### **QCD fit function**:
# 
# $$ \frac{d\sigma}{dp_T} = A \cdot \left(1 + \frac{p_T}{p_0} \right)^{-n} $$
# 

# 
