"""Rebuild the presentation figures from the three supplied files.
Python dependencies: numpy, pandas, openpyxl (read only), matplotlib.
Run: python make_tls_figures.py
All calculations are descriptive operational summaries, not validation metrics.
"""
from pathlib import Path
from io import BytesIO
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
w=pd.read_excel(ROOT/'inputs/StemSegm_CIrcumference.xlsx').sort_values('TreeID')
s=pd.read_csv(ROOT/'inputs/trunk_trace_summary.csv').sort_values('TreeID')
d=pd.read_csv(ROOT/'inputs/trunk_slice_checks.csv').sort_values(['TreeID','Bin'])
for name in ['Accepted','PathUsesCoverage']:
 if d[name].dtype!=bool:
  d[name]=d[name].astype(str).str.lower().map({'true':True,'false':False})
  assert not d[name].isna().any(),name
assert not d.duplicated(['TreeID','Bin']).any()
assert set(s.TreeID)==set(w.TreeID)==set(d.TreeID)
a=d[d.Bin==0].sort_values('TreeID')
assert np.allclose(w.DBH,a.Radius_m*200)
assert np.allclose(s.PilotDiameterAt1p3_cm,w.DBH)
assert np.allclose(w['Stem Circumference'],w.DBH*3.14159)
for _,r in s.iterrows():
 g=d[d.TreeID==r.TreeID]
 for col,n in [('AcceptedSlices',g.Accepted.sum()),('FailedSlices',(~g.Accepted).sum()),('StrictSlices',g.FitQuality.eq(1).sum()),('CoverageReviewSlices',g.FitQuality.eq(2).sum()),('PathUsesCoverageSlices',(g.Accepted & g.PathUsesCoverage).sum())]:
  assert r[col]==n,(r.TreeID,col)
assert (d.Accepted==d.FitQuality.isin([1,2])).all()
# Four disjoint categories preserve the distinction between local fit quality
# and a trace whose route already used a provisional coverage slice.
d['PlotClass']=np.select([~d.Accepted,d.FitQuality.eq(2),d.PathUsesCoverage],['Rejected','Coverage trial','Strict after coverage'],'Strict before coverage')
acc=d[d.Accepted].copy();bad=d[~d.Accepted].copy()
NAVY='#15354B'; TEAL='#167D9A'; LIGHT='#8CC7D4'; GOLD='#D18A14'; RED='#BC4C4C'; GREY='#667783'; GRID='#E4EBEF'
C={'Strict before coverage':TEAL,'Strict after coverage':LIGHT,'Coverage trial':GOLD,'Rejected':RED}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':15,'axes.labelsize':16,'axes.titlesize':18,'axes.titleweight':'bold','text.color':NAVY,'axes.labelcolor':NAVY,'xtick.color':NAVY,'ytick.color':NAVY,'axes.edgecolor':'#A9B7C1','axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','axes.facecolor':'white','pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
pdf=PdfPages(ROOT/'TLS_Operational_Plots.pdf')
figures=[]
def canvas(title,sub):
 f=plt.figure(figsize=(16,9))
 f.text(.07,.932,title,fontsize=27,fontweight='bold',ha='left')
 f.text(.07,.875,sub,fontsize=15,color=GREY,ha='left')
 return f

def foot(f,text):f.text(.07,.045,text,fontsize=11.5,color=GREY,va='bottom')
def clean(ax,axis='x'):
 ax.set_axisbelow(True);ax.grid(axis=axis,color=GRID,lw=.9)
 ax.tick_params(length=0,pad=8)
 ax.spines['left'].set_visible(False)

def save(f,name):
 buf=BytesIO(); f.savefig(buf,format='png',dpi=240,facecolor='white')
 (OUT/(name+'.png')).write_bytes(buf.getvalue())
 f.savefig(OUT/(name+'.svg'),facecolor='white')
 pdf.savefig(f,facecolor='white')
 figures.append((name,f))

# 01 Stem size. Circumference is derived, never represented as independent data.
f=canvas('Fitted stem diameters span 21.1–44.4 cm','Cross-sections at 1.3 m above ground • seven trees in the 10 × 10 m pilot')
ax=f.add_axes([.11,.19,.61,.57]);y=np.arange(len(w))
ax.barh(y,w.DBH,color=TEAL,height=.57)
for yy,diam,circ in zip(y,w.DBH,w['Stem Circumference']):
 ax.text(diam+.6,yy,f'{diam:.1f}',va='center',fontweight='bold',fontsize=16)
 ax.text(1.21,yy,f'{circ:.1f}',transform=ax.get_yaxis_transform(),ha='center',va='center',fontsize=16)
ax.set(yticks=y,yticklabels=[f'Tree {i}' for i in w.TreeID],xlim=(0,52),xlabel='Fitted diameter at 1.3 m (cm)')
ax.invert_yaxis();clean(ax)
f.text(.847,.80,'Derived circumference',fontsize=16,fontweight='bold',ha='center')
f.text(.847,.765,'cm',fontsize=14,color=GREY,ha='center')
foot(f,'Circumference = π × fitted diameter; the workbook uses π = 3.14159.\nThese are horizontal circle estimates from the sampled TLS pilot, not independent field measurements.')
save(f,'01_Stem_Diameter_and_Circumference')

# 02 Positions from the fitted 1.3 m cross-sections. No CRS is assigned.
f=canvas('Seven fitted stem centres within the pilot','Positions at 1.3 m above ground • circle area is proportional to fitted diameter²')
ax=f.add_axes([.10,.18,.64,.62]);origin=np.array([420402.6,6481810.0])
x=w.X.to_numpy()-origin[0]; yy=w.Y.to_numpy()-origin[1];sizes=(w.DBH/30)**2*550
ax.scatter(x,yy,s=sizes,c=TEAL,alpha=.85,edgecolor='white',linewidth=2,zorder=3)
for xx,yyy,tid,diam in zip(x,yy,w.TreeID,w.DBH):
 ax.annotate(f'Tree {tid}\n{diam:.1f} cm',(xx,yyy),xytext=(12,-6 if yyy>9 else 10),textcoords='offset points',fontsize=13,fontweight='bold',va='top' if yyy>9 else 'bottom')
ax.set(xlim=(0,10),ylim=(0,10),aspect='equal',xticks=np.arange(0,11,2),yticks=np.arange(0,11,2),xlabel='Local X (m)',ylabel='Local Y (m)');clean(ax,'both')
for diameter in [20,30,40]:ax.scatter([],[],s=(diameter/30)**2*550,c=TEAL,label=f'{diameter} cm')
ax.legend(title='Fitted diameter',bbox_to_anchor=(1.22,.62),loc='center left',frameon=False,labelspacing=1.8)
f.text(.76,.22,'Pilot boundary\n10 × 10 m',fontsize=18,fontweight='bold')
foot(f,'Local origin: X = 420402.6 m, Y = 6481810.0 m in the input coordinates; CRS unspecified.\nSymbols are enlarged for visibility and do not show actual trunk footprints or crowns.')
save(f,'02_Pilot_Stem_Locations')

# 03 Retained height and point counts, neither interpreted as total tree size.
f=canvas('Retained trunk points reach 13.4–18.2 m above ground','85,283 candidate trunk points across seven trees • retained extent, not total tree height')
left=f.add_axes([.11,.21,.40,.53]);right=f.add_axes([.65,.21,.28,.53]);y=np.arange(7)
for ax in [left,right]:ax.set(yticks=y,yticklabels=[f'Tree {i}' for i in s.TreeID]);ax.invert_yaxis();clean(ax)
left.hlines(y,s.TracedMinHAG_m,s.TracedMaxHAG_m,color=LIGHT,lw=9)
left.scatter(s.TracedMaxHAG_m,y,s=125,color=TEAL,zorder=4)
for yy,h in zip(y,s.TracedMaxHAG_m):left.text(h+.4,yy,f'{h:.1f} m',va='center',fontweight='bold',fontsize=14)
left.set(xlim=(0,21),xlabel='Highest retained point (m above ground)');left.set_title('Vertical extraction extent',loc='left',pad=22)
right.barh(y,s.CandidatePoints/1000,color=TEAL,height=.55)
for yy,n in zip(y,s.CandidatePoints):right.text(n/1000+.25,yy,f'{n:,}',va='center',fontsize=13,fontweight='bold')
right.set(xlim=(0,22),xticks=[0,5,10,15,20],xlabel='Candidate points (thousands)',yticklabels=[]);right.set_title('Retained point count',loc='left',pad=22)
foot(f,'Every tree starts at approximately 0.20 m above ground. All upward traces stop after two consecutive failed slices.\nA high point count is sampling support, not a measure of extraction accuracy or completeness.')
save(f,'03_Trunk_Extent_and_Point_Count')

# 04 Exact slice statuses with no interpolation into untested heights.
f=canvas('540 of 563 tested slices were retained','528 passed the original checks • 12 used the coverage trial • 23 were rejected')
ax=f.add_axes([.10,.25,.84,.56])
for _,r in d.iterrows():ax.add_patch(Rectangle((r.TreeID-.34,r.Height_m-.10),.68,.20,facecolor=C[r.PlotClass],edgecolor='white',linewidth=.15))
ax.set(xlim=(.4,7.6),ylim=(0,20),xticks=range(1,8),xticklabels=[f'Tree {i}' for i in range(1,8)],yticks=range(0,21,2),ylabel='Slice centre height above ground (m)');clean(ax,'y')
legend=[Patch(facecolor=C[k],label=label) for k,label in [('Strict before coverage','Strict, no coverage dependency (429)'),('Strict after coverage','Strict, after a coverage trial (99)'),('Coverage trial','Coverage trial accepted (12)'),('Rejected','Rejected (23)')]]
f.legend(handles=legend,loc='lower center',bbox_to_anchor=(.53,.11),ncol=2,frameon=False,fontsize=13)
foot(f,'Each band represents a tested 0.2 m slice. Blank space above a column was not tested.\nThe 99 pale-blue slices pass local strict checks, but their tracing path previously used a provisional coverage slice.')
save(f,'04_Slice_Status_by_Height')

# 05 Diameter profiles. Lines break at rejected slices; rejected radii omitted.
f=canvas('Accepted circle fits describe the extracted trunk profiles','Horizontal fitted diameters along the retained stem sections • shared axes across all seven trees')
axes=f.subplots(2,4,sharex=True,sharey=True);f.subplots_adjust(left=.075,right=.96,bottom=.15,top=.76,hspace=.38,wspace=.23)
for ax,tid in zip(axes.flat,range(1,8)):
 g=acc[acc.TreeID==tid].sort_values('Bin')
 for _,run in g.groupby(g.Bin.diff().ne(1).cumsum()):ax.plot(200*run.Radius_m,run.Height_m,color='#CCD6DD',lw=1.7,zorder=1)
 for cls in C:
  z=g[g.PlotClass==cls]
  if len(z):ax.scatter(200*z.Radius_m,z.Height_m,s=17 if cls!='Coverage trial' else 42,color=C[cls],edgecolor='white' if cls=='Coverage trial' else 'none',linewidth=.5,zorder=3)
 ax.axhline(1.3,color=GREY,ls=':',lw=1)
 ax.set(xlim=(0,55),ylim=(0,20),xticks=[0,20,40],yticks=[0,5,10,15,20]);clean(ax,'both')
 ax.set_title(f'Tree {tid}',loc='left',fontsize=17,pad=9)
axes.flat[7].set_axis_off();key=axes.flat[7]
key.text(.03,.92,'Reading the profiles',transform=key.transAxes,fontsize=16,fontweight='bold')
for pos,cls,label in [(.74,'Strict before coverage','Strict fit'),(.58,'Strict after coverage','Strict after coverage'),(.42,'Coverage trial','Coverage trial')]:
 key.plot(.07,pos,marker='o',markersize=6,color=C[cls],transform=key.transAxes)
 key.text(.15,pos,label,transform=key.transAxes,fontsize=12,va='center')
key.text(.03,.015,'Dotted line = 1.3 m\nGaps are not interpolated.\nRejected fits are omitted.',transform=key.transAxes,fontsize=10.5,color=GREY,linespacing=1.15,va='bottom')
f.text(.50,.098,'Fitted diameter (cm)',ha='center',fontsize=16);f.text(.018,.48,'Height above ground (m)',rotation=90,va='center',fontsize=16)
foot(f,'Profiles are unsmoothed horizontal circle fits. Local widening can reflect stem shape, lean, branch junctions or fitting effects.')
save(f,'05_Trunk_Diameter_Profiles')

# 06 Internal fitting diagnostics; no comparison with independent truth.
f=canvas('Retained slices have a median radial residual of 5.3 mm','Internal fitting diagnostics for 540 accepted slices • points show individual slices, boxes show their distributions')
ax1=f.add_axes([.10,.23,.36,.53]);ax2=f.add_axes([.59,.23,.36,.53]);rng=np.random.default_rng(20260925)
for ax,col,scale,limit,threshold in [(ax1,'RadialRMSE_m',1000,(0,14),12),(ax2,'RingFraction',100,(0,105),35)]:
 vals=[acc.loc[acc.TreeID==tid,col].to_numpy()*scale for tid in range(1,8)]
 bp=ax.boxplot(vals,positions=np.arange(1,8),widths=.48,showfliers=False,patch_artist=True,medianprops={'color':NAVY,'linewidth':2},boxprops={'facecolor':'#DCECF0','edgecolor':TEAL},whiskerprops={'color':GREY},capprops={'color':GREY})
 for tid in range(1,8):
  g=acc[acc.TreeID==tid];jit=rng.uniform(-.13,.13,len(g));cols=np.where(g.FitQuality.eq(2),GOLD,TEAL)
  ax.scatter(tid+jit,g[col]*scale,s=13,c=cols,alpha=.57,zorder=3)
 ax.axhline(threshold,color=RED,ls='--',lw=1.5)
 ax.set(xlim=(.4,7.6),ylim=limit,xticks=range(1,8),xlabel='Tree ID');clean(ax,'y')
ax1.set(ylabel='Radial RMSE of ring inliers (mm)');ax1.set_title('Circle-fit residual',loc='left',pad=22)
ax1.text(.02,.92,'Original limit 12 mm',transform=ax1.transAxes,color=RED,fontsize=12)
ax2.set(ylabel='Ring inliers / search-window points (%)');ax2.set_title('Support within the search window',loc='left',pad=22)
ax2.text(.02,.36,'Minimum fraction 35%',transform=ax2.transAxes,color=RED,fontsize=12)
f.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=TEAL,label='Original checks passed'),Line2D([0],[0],marker='o',color='none',markerfacecolor=GOLD,label='Coverage trial')],loc='lower center',bbox_to_anchor=(.52,.11),ncol=2,frameon=False,fontsize=13)
foot(f,'RMSE measures distances to the fitted circle; it is not independent DBH accuracy. Coverage trials use the tighter 8 mm RMSE limit.\nBox = interquartile range; line = median; whiskers extend to observations within 1.5 × IQR. Only retained slices are shown.')
save(f,'06_Internal_Fit_and_Point_Support')

# 07 Multi-response rejection diagnostics from the previous-circle model.
counts=bad.StrictChecksFailed.fillna('').str.split(';').explode();counts=counts[counts.ne('')].value_counts()
labels={'angular_gap_limit':'Angular gap above 120°','ring_points_below_min':'Fewer than 20 ring inliers','occupied_sectors_below_min':'Fewer than 8 occupied sectors','ring_fraction_below_min':'Ring fraction below 35%','search_points_below_min':'Fewer than 20 search points','prior_fit_optim_error':'Optimizer error'}
f=canvas('Angular gaps and limited ring support dominate the rejections','Recorded conditions for the previous-circle model in 23 rejected slices • multiple conditions can occur in one slice')
ax=f.add_axes([.30,.21,.57,.55]);y=np.arange(len(counts))
ax.barh(y,counts.to_numpy(),color=[RED]*3+[GREY]*(len(counts)-3),height=.58)
for yy,n in zip(y,counts):ax.text(n+.3,yy,f'{n} '+('slice' if n==1 else 'slices'),va='center',fontsize=16,fontweight='bold')
ax.set(yticks=y,yticklabels=[labels.get(k,k) for k in counts.index],xlim=(0,18),xticks=range(0,19,3),xlabel='Rejected slices with the recorded condition');ax.invert_yaxis();clean(ax)
foot(f,'There are 9 internal rejected slices and 14 terminal rejections (two per tree). Counts overlap and must not be summed as distinct slices.\nOne slice also logged optimizer code 52 (Tree 3, 11.1 m). These logs do not report the reason for every RANSAC proposal.')
save(f,'07_Rejected_Slice_Conditions')
pdf.close()
# Data tables and concise presentation notes for traceability.
summary=s.merge(w[['TreeID','Stem Circumference']],on='TreeID',validate='one_to_one')
summary.rename(columns={'Stem Circumference':'DerivedCircumference_cm'},inplace=True)
summary['TestedSlices']=summary.AcceptedSlices+summary.FailedSlices
summary['RetainedFraction_tested']=summary.AcceptedSlices/summary.TestedSlices
summary['StrictAfterCoverageSlices']=summary.PathUsesCoverageSlices-summary.CoverageReviewSlices
summary.to_csv(ROOT/'Operational_Summary.csv',index=False)
notes='''# TLS pilot operational figures

These figures describe the seven-tree central 10 × 10 m pilot. They do not constitute independent validation, a full-plot inventory, branch segmentation, or total tree-height estimation.

## What the three files contain

- `StemSegm_CIrcumference.xlsx` contains the seven fitted cross-sections at 1.3 m. Its `DBH` values are in centimetres, calculated as 200 × radius in metres. `Stem Circumference` is also in centimetres and is computed using 3.14159 × DBH. It is not a separate measurement.
- `trunk_trace_summary(2).csv` contains per-tree retained point counts, accepted/rejected slice counts, retained vertical extent and stopping conditions.
- `trunk_slice_checks(1).csv` contains all 563 tested slices, with local fit status, supporting points, residuals, angular coverage and coverage-path flags.

All seven diameter values and all per-tree slice counts agree between these inputs. The circumference calculation was retained as supplied; replacing 3.14159 with full precision pi would change the largest value by less than 0.00012 cm.

## Operational findings

- Seven trees and 85,283 candidate trunk points.
- 540 retained slices out of 563 tested (95.9%). This is a retention rate among attempted slices, not segmentation accuracy or whole-stem completeness.
- 528 strict fits, 12 provisional coverage-trial fits, and 23 rejections.
- Of the 528 strict fits, 99 occur on a route that previously used a coverage trial. Including the 12 trial slices gives 111 retained slices flagged `PathUsesCoverage`.
- The fitted diameter at 1.3 m ranges from 21.1 cm (Tree 2) to 44.4 cm (Tree 5).
- The highest retained point ranges from 13.399 m (Tree 3) to 18.189 m (Tree 7). Tree 4 is nearly equal at 18.183 m. This is retained trunk extent, not total tree or trunk height.
- All seven traces end after two consecutive rejected slices. There are nine earlier, internal rejections.
- Median accepted-slice radial RMSE is 5.32 mm. This is an internal circle-fit residual, not independent diameter accuracy.
- At 1.3 m, all seven fits pass the original checks and occupy all 12 angular sectors.
- Rejected model conditions include angular gaps in 15 slices, fewer than 20 ring points in 14, and fewer than 8 occupied sectors in 13. Conditions overlap. One rejected slice (Tree 3 at 11.1 m) logs optimizer code 52; no fitting-error cause is inferred beyond that log.

## Suggested speaking notes by figure

1. **Stem diameter and circumference.** Tree 5 has the largest fitted stem cross-section and Tree 2 the smallest. Circumference is a conversion of the same fitted diameter.
2. **Pilot stem locations.** These are fitted stem centres at 1.3 m, plotted within the pilot boundary. Enlarged symbols communicate diameter, not trunk footprints. No CRS is inferred.
3. **Trunk extent and point count.** Trunk extraction extends above 13 m for every tree; Tree 4 has the most retained points. Point density and retained height describe different aspects of the output.
4. **Slice status by height.** Most tested slices meet the original criteria. Amber bands mark limited coverage trials; pale blue bands show strict fits whose route depends on earlier trials. White space above each column was not tested.
5. **Trunk diameter profiles.** Accepted horizontal-circle diameters generally narrow upwards, with local changes. Gaps are left open and rejected radii are not plotted as valid stem sizes. These profiles have not been smoothed or corrected for lean.
6. **Internal fit and point support.** This shows how closely the retained ring points fit circles and what fraction of searched points supports each circle. Neither diagnostic is a validation score.
7. **Rejected slice conditions.** Angular gaps and sparse ring support are common operational constraints. These are logs of the previous-circle candidate; rejected RANSAC proposals are not individually diagnosed.

## Use and reproduce

Each figure is 16:9, with a 3840 × 2160 pixel PNG and a vector SVG. The combined PDF contains one figure per page. SVG text remains text, subject to font availability in the destination application.

Run `python make_tls_figures.py` from this extracted folder to regenerate the figures, PDF and summary. Requires Python with numpy, pandas, openpyxl (for reading the source workbook) and matplotlib. Source inputs are included unchanged. The script uses a fixed seed only to position the diagnostic dots horizontally; no source measurements are randomized.
'''
(ROOT/'README_Presentation_Notes.md').write_text(notes)
# A lightweight contact sheet for visual inspection, not a presentation figure.
from PIL import Image,ImageOps,ImageDraw
thumbs=[]
for name,_ in figures:
 im=Image.open(OUT/(name+'.png')).convert('RGB');im.thumbnail((800,450))
 tile=Image.new('RGB',(800,485),'#eef2f4');tile.paste(im,(0,0));ImageDraw.Draw(tile).text((10,460),name,fill='#15354B');thumbs.append(tile)
sheet=Image.new('RGB',(1600,4*485),'white')
for i,tile in enumerate(thumbs):sheet.paste(tile,((i%2)*800,(i//2)*485))
sheet.save(ROOT/'contact_sheet.jpg',quality=88)
print(summary[['TreeID','CandidatePoints','PilotDiameterAt1p3_cm','DerivedCircumference_cm','TracedMaxHAG_m','AcceptedSlices','FailedSlices','CoverageReviewSlices','StrictAfterCoverageSlices']].to_string(index=False))
print('Wrote',len(figures),'figures (PNG/SVG), combined PDF, summary and notes.')
plt.close('all')
