'''
此文件用来绘制绘制实验结果图
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, auc
from sklearn.utils import resample
from math import pi
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['svg.fonttype'] = 'none'

# --------------Panel A ROC_curve------------#
def plot_roc_curve(data, model_name, color, alpha=1.0):
    data = pd.read_csv(data)
    all_fpr, all_tpr = data['fpr'].values, data['tpr'].values
    all_auc = auc(all_fpr, all_tpr)
    plt.plot(all_fpr, all_tpr, color=color, lw=2, alpha=alpha,
            label=f'{model_name} AUC = {all_auc:.3f}')
    return

plt.figure(figsize=(4, 3.6))
ours_path = "table/Ours_overall_ROC_data.csv"
radiomics_path = "table/Radiomics_overall_ROC_data.csv"
resnet18_path = 'table/Resnet18_overall_ROC_data.csv'
plot_roc_curve(resnet18_path, 'ResNet-18', '#1F77B4', 0.5)
plot_roc_curve(radiomics_path, 'Radiomics', '#ffc477', 0.8)
plot_roc_curve(ours_path, 'CT4CMSPlus', '#f08877', 1.0)

plt.plot([0, 1], [0, 1], '--', color='#AFABAB', lw=1)
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('1 - Specificity', fontsize=14, labelpad=10)
plt.ylabel('Sensitivity', fontsize=14, labelpad=10)
plt.legend(loc="lower right", fontsize=11)
plt.tight_layout()
plt.savefig(f'img_plot/ROC_curve.png', dpi=300)
plt.savefig("img_plot/ROC_curve.svg", format="svg", bbox_inches="tight")
plt.savefig("img_plot/ROC_curve.pdf", format="pdf", bbox_inches="tight", dpi=300)
plt.show()

#--------------Panel B confusion_matrix-----------#
df = pd.read_csv("table/Ours_five_fold_results.csv")
cm = confusion_matrix(df['labels'], df['preds'])
classes = ["CMS4-TME+", "CMS4-TME-"]

accuracy = np.trace(cm) / np.sum(cm)  
print(f"Accuracy: {accuracy:.3f}")

plt.figure(figsize=(3, 3))
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['svg.fonttype'] = 'none'
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", cbar=True,
            xticklabels=classes, yticklabels=classes, annot_kws={"size": 12}, cbar_kws={"shrink": 0.8}, vmin=0, vmax=cm.max() * 1.2, square=True)

# plt.title("Confusion Matrix", fontsize=12)
plt.xlabel("Prediction", fontsize=10)
plt.ylabel("Reference", fontsize=10)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_aspect('equal')  
plt.text(0.5, -0.32, f"Accuracy: {accuracy:.3f}", fontsize=10, ha="center", transform=plt.gca().transAxes)

plt.tight_layout()
plt.savefig(f'img_plot/confusion_matrix.png', dpi=300)
plt.savefig("img_plot/confusion_matrix.svg", format="svg", bbox_inches="tight")
plt.savefig("img_plot/confusion_matrix.pdf", format="pdf", bbox_inches="tight")
plt.show()

#--------------Panel C Radar Chart------------#
def metrics(file, auc):

    all_true = file['labels'].apply(lambda x: 0 if x == 'CMS0' else 1).values
    all_pred = file['preds'].apply(lambda x: 0 if x == 'CMS0' else 1).values
    ours_auc = auc

    tn, fp, fn, tp = confusion_matrix(all_true, all_pred).ravel()
    acc = (tp + tn) / (tp + tn + fp + fn) 
    sen = tp / (tp + fn) if (tp + fn) != 0 else 0  
    spe = tn / (tn + fp) if (tn + fp) != 0 else 0  
    ppv = tp / (tp + fp) if (tp + fp) != 0 else 0  
    npv = tn / (tn + fn) if (tn + fn) != 0 else 0  
    print(f"Accuracy (ACC): {acc:.3f}")
    print(f"Sensitivity (SEN): {sen:.3f}")
    print(f"Specificity (SPE): {spe:.3f}")
    print(f"Positive Predictive Value (PPV): {ppv:.3f}")
    print(f"Negative Predictive Value (NPV): {npv:.3f}")
    return [auc, acc, sen, spe, ppv, npv]

categories = ["AUC", "ACC", "SEN", "SPE", "PPV", 'NPV']
model_names = ["Ours", "Radiomics", "ResNet-18"]

ours_path = pd.read_csv("table/Ours_five_fold_results.csv")
radiomics_path = pd.read_csv("table/Radiomics_five_fold_results.csv")
resnet18_path = pd.read_csv("table/Resnet18_five_fold_results.csv")

ours_auc = 0.887  
radiomics_auc = 0.802
resnet18_auc = 0.794
Ours = metrics(ours_path, ours_auc)
Radiomics = metrics(radiomics_path, radiomics_auc)
ResNet18 = metrics(resnet18_path, resnet18_auc)

N = len(categories)
angles = [n / float(N) * 2 * pi for n in range(N)]

angles += angles[:1]
categories += categories[:1] 

Radiomics += Radiomics[:1]
ResNet18 += ResNet18[:1]
Ours += Ours[:1]

fig, ax = plt.subplots(figsize=(5, 4), subplot_kw={"polar": True})
ax.set_theta_offset(pi / 2)
ax.set_theta_direction(-1)

ax.fill(angles, Radiomics, color="#ffc477", alpha=0.1)
ax.plot(angles, Radiomics, color="#ffc477", linewidth=1.5, alpha=0.9, label="Radiomics")

ax.fill(angles, ResNet18, color="#1F77B4", alpha=0.1)
ax.plot(angles, ResNet18, color="#1F77B4", linewidth=1.5, alpha=0.5,label="ResNet-18")

ax.fill(angles, Ours, color="#f08877", alpha=0.1)
ax.plot(angles, Ours, color="#f08877", linewidth=1.5, alpha=1,label="CT4CMSPlus")

ax.set_rlabel_position(180)
ax.set_yticks([0.65, 0.75, 0.85, 0.9]) 
ax.set_yticklabels(["0.65", "0.75", "0.85", "0.90"], fontsize=8, color="#3C3C3C")#fontweight="bold"
ax.set_ylim(0.65, 0.9)
ax.grid(color="gray", linestyle="--", linewidth=0.5)
ax.set_xticks(angles)  
ax.set_xticklabels(categories, fontsize=12)
ax.legend(
    loc="upper center",  
    bbox_to_anchor=(0.5, 1.3),  
    fontsize=11, 
    ncol=3,  
    frameon=False,
    # title="Models",  
    title_fontsize=12  
)

plt.rcParams['font.family'] = 'Arial'
plt.tight_layout()
plt.savefig(f'img_plot/Radar.png', dpi=300)
plt.savefig("img_plot/Radar.svg", format="svg", bbox_inches="tight")
plt.savefig("img_plot/Radar.pdf", format="pdf", bbox_inches="tight", dpi=300)
plt.show()

#--------------Panel E Immunity Response Barplot------------#
data = {
    'Response': ['iPR', 'iSD', 'iUPD', 'iCPD'],
    'CMS1': [3, 7, 1, 0],
    'CMS2': [3, 7, 4, 4],
    'CMS3': [5, 9, 1, 5],
    'CMS4-TME+': [10, 7, 1, 4],
    'CMS4-TME-': [7, 13, 1, 5],
}
df = pd.DataFrame(data)
df['Freq_CMS1'] = df['CMS1'] / df['CMS1'].sum()
df['Freq_CMS2'] = df['CMS2'] / df['CMS2'].sum()
df['Freq_CMS3'] = df['CMS3'] / df['CMS3'].sum()
df['Freq_CMS4-TME+'] = df['CMS4-TME+'] / df['CMS4-TME+'].sum()
df['Freq_CMS4-TME-'] = df['CMS4-TME-'] / df['CMS4-TME-'].sum()

fig, ax = plt.subplots(figsize=(4, 4))
bar_width = 0.15
index = list(range(len(df['Response']))) 

bar1 = ax.bar(index, df['Freq_CMS1'], bar_width, label='CMS1', color='#E69F00', alpha=0.6)
bar2 = ax.bar([p + bar_width for p in index], df['Freq_CMS2'], bar_width, label='CMS2', color="#377EB8", alpha=0.6)
bar3 = ax.bar([p + bar_width * 2 for p in index], df['Freq_CMS3'], bar_width, label='CMS3', color="#CC7AA7", alpha=0.6)
bar4 = ax.bar([p + bar_width * 3 for p in index], df['Freq_CMS4-TME+'], bar_width, label='CMS4-TME+', color='#3D5588', alpha=0.6)
bar5 = ax.bar([p + bar_width * 4 for p in index], df['Freq_CMS4-TME-'], bar_width, label='CMS4-TME-', color="#8491B5", alpha=0.6)

# 设置图形标题和标签
ax.set_xlabel('Response category', fontsize=14)
ax.set_ylabel('Response (rel. freq.)', fontsize=14)
#ax.set_title('Response to Cetuximab Monotherapy', fontsize=14)
ax.set_xticks([p + bar_width * 1.5 for p in index])  # 调整 X 轴刻度位置
ax.set_xticklabels(df['Response'])
ax.legend(fontsize=10, loc='upper right')
# 移除上框线和右框线
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# 显示图形
plt.tight_layout()
# plt.savefig(f'img_plot/immunity_freq.png', dpi=300)
plt.savefig("img_plot/immunity_freq.svg", format="svg", bbox_inches="tight")
plt.savefig("img_plot/immunity_freq.pdf", format="pdf", bbox_inches="tight", dpi=300)
plt.show()

#--------------Panel F Immunity Response Barplot------------#
categories = ['CMS1/4-TME+', 'CMS2/3/4-TME-']
no_responder = [6/(6+28), 20/(35+20)]  
responder = [28/(28+6), 35/(20+35)]    

bar_width = 0.1
x = np.arange(len(categories)) * 0.12  

plt.figure(figsize=(2, 3.6)) 
plt.bar(x, responder, color="#4F9292", width=bar_width, label='Responder')  
plt.bar(x, no_responder, color="#D8A76F", width=bar_width, bottom=responder, label='Non-responder') 
plt.text(np.mean(x), 1.1, r'$\mathit{P} = 4.78e{-}02$', ha='center', fontsize=10)

plt.xticks(x, categories, rotation=30, fontsize=10)  
plt.ylabel('Percentage', fontsize=12) 
plt.ylim(0, 1.05)
plt.yticks(np.arange(0, 1.1, 0.25))  
#plt.legend(title='Response', loc='upper center', fontsize=10) 

plt.tight_layout()
plt.savefig(f'img_plot/immunity_plot.png', dpi=300)
plt.savefig("img_plot/immunity_plot.svg", format="svg", bbox_inches="tight")
plt.savefig("img_plot/immunity_plot.pdf", format="pdf", bbox_inches="tight", dpi=300)
plt.show()
