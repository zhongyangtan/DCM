import json
import os
import codecs
import re
import pandas as pd
import openpyxl
import chardet  
import time
from joblib import Parallel, delayed
import multiprocessing
import numpy as np



ex_num = 100
ex_den = 600

flag_num = 20  
map_add = 0   

count=0
count_txt = 0
count_csv = 0
count_mer = 0
txt_c = 0 

flag_name = 0  
peak_name_new = '' 



    
with open("config_segname.json", "r") as fp_seg:
    config_seg = json.loads(fp_seg.read())

with open("config_peak.json", "r") as jn:
    config_peak = json.loads(jn.read())
    
with open("config_control.json", "r") as jn_control:
    config_control = json.loads(jn_control.read())
 
# type_1 = 1
fileInPath = './ssr_in/'
#fileOutPath = './ssr_density_statistics/'
ssr_allPath = './ssr_all/'
fileHead = 'Type,position,label,ratio_mono,ratio_di,ratio_tri,ratio_tetra,' \
               'ratio_penta,ratio_hexa,ratio_total,ratio_len_mono,ratio_len_di,ratio_len_tri,' \
               'ratio_len_tetra,ratio_len_penta,ratio_len_hexa,ratio_len_total,ex_num,ex_den\n'

#path = os.getcwd().replace('\\', '/')

folders = ['map_den','ssr_in', 'ssr_all','infile','xlsx_to_json','in_outfile','lastfile','lastfile/LP','ssr_all/ALL_mer','ucsc_bed','ucsc_bed/bedGf','ucsc_bed/region']

for folder in folders:
    if not os.path.isdir(folder):
        os.mkdir(folder)

print("\n")
print("Processing, please wait 》》》》》")
print("\n")

infile = r"infile\\"
outfile =  r"ssr_in\\" 
all_file = r"ssr_all\\"
in_outfile = r"in_outfile\\"
xlsx_to_json =r"xlsx_to_json\\"
temp_map = r"map_den\\"
lp_save = "lastfile/LP"
ALL_mer = 'ALL_mer\\'
mer_file = 'ssr_all\\ALL_mer\\'
bed_outfile = "ucsc_bed\\"
land_bgf = 'ucsc_bed\\bedGf\\'
seg_file ="seg_file\\"


flag_replace = '1'

seg_list = os.listdir(seg_file)  

global data_seg

if len(seg_list) != 0:
    data_seg = pd.read_excel(seg_file + seg_list[0], engine="openpyxl")  
    old_name = 'G38'
    new_name = 'GRCh38'
else:
    old_name = 'HS'
    new_name = 'CHM13'


deal_file = os.listdir(xlsx_to_json) 
df = pd.read_excel(xlsx_to_json+deal_file[0],engine="openpyxl") 
df1 = df.drop(['id', 'to'], axis=1) 
df1['from'] = df1['from'].apply(lambda x: x - 1) 
df1.rename(columns={'from': 'num'}, inplace=True) 
df1['num'] = df1['num'].astype(int)  
df1.to_json("config_txt.json",orient='records') 

time.sleep(1)

with open("config_txt.json", "r") as fp_txt:
    config_txt = json.loads(fp_txt.read())


df2 = df.drop(['id'], axis=1)
df2.rename(columns={'from': 'start','to': 'end','filename': 'name'}, inplace=True)
name_seg = config_seg[0]["segname1"]
#print(name_seg)
df2['segname'] = df2.apply(lambda x: name_seg+"_chr"+x['name'].split('-')[1]+": "+x['name'].split('-S')[0], axis=1)   
df2 = df2.reindex(columns=['name','start','end','segname'])  
df2['start'] = df2['start'].astype(int)
df2['end'] = df2['end'].astype(int)
df2.to_excel(temp_map+deal_file[0].split('.')[0]+"-json"+".xlsx",index=False)
df2.to_json("config.json",orient='records')

time.sleep(1)
with open("config.json", "r") as fp:
    config = json.loads(fp.read())
    
def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def get_min_repr(seq):
    sort_len = len(seq)
    if sort_len == 1:
        return seq
    elif sort_len in [2, 3, 4, 5, 6]:
        seq_numbers = seq.replace('A', '1').replace('C', '2').replace('G', '3').replace('T', '4')
    else:
        raise ValueError("Data Processing Error---check_motif")

    sort_all = seq_numbers + seq_numbers
    sort_all_len = len(sort_all)
    last_sort = []
    count = sort_len

    for i in range(sort_len):
        last_sort.append(sort_all[i:(sort_all_len - count)]) 
        count -= 1

    sorted_repr = sorted(last_sort, key=lambda i: int(re.match(r'(\d+)', i).group()))
    min_repr = sorted_repr[0].replace('1', 'A').replace('2', 'C').replace('3', 'G').replace('4', 'T')
    sort_hexa = "TTAGGGTTAGGG"
    if min_repr in sort_hexa:
        min_repr ="TTAGGG"
    #return sort_temp
    return min_repr

results_df = pd.DataFrame(columns=['Orig_Seq', 'Adj_Seq', 'Comp_Pal_Seq'])

def process_sequence_rev(seq): 
    if seq in "TTAGGGTTAGGG" and len(seq) == 6:
        return seq
    if seq in "AACCCTAACCCT" and len(seq) == 6:
        seq_deal = 'CCCTAA'
        min_rev_seq = 'TTAGGG'
        #min_seq = reverse_complement(seq)
        if seq not in results_df['Orig_Seq'].values:
            results_df.loc[len(results_df)] = [seq, seq_deal, min_rev_seq]
        return seq_deal
    rev_seq = reverse_complement(seq)
    #print(rev_seq)
    min_seq = get_min_repr(rev_seq)
    #print(min_seq)
    #min_rev_seq = get_min_repr(rev_seq)

    seq_numbers = seq.replace('A', '1').replace('C', '2').replace('G', '3').replace('T', '4')
    min_seq_numbers = min_seq.replace('A', '1').replace('C', '2').replace('G', '3').replace('T', '4')

    if seq_numbers <= min_seq_numbers:
        return seq
    else:
        min_rev_seq = reverse_complement(min_seq)
        rev_min = reverse_complement(min_rev_seq)
        if rev_min != min_seq:
            print("Error---seq_rev")
        if seq == min_rev_seq:
            return seq
        else:
            if seq not in results_df['Orig_Seq'].values:
                results_df.loc[len(results_df)] = [seq, min_rev_seq, min_seq]
            return  min_rev_seq
    
    
def check_motif(motif_param): 
    last_sort = []
    sort_len = len(motif_param)
    if sort_len == 2 or sort_len == 4 or sort_len ==5 or sort_len==6 or sort_len == 3: 
        sort_temp = motif_param.replace('A', '1')
        sort_temp = sort_temp.replace('C', '2')
        sort_temp = sort_temp.replace('G', '3')
        sort_temp = sort_temp.replace('T', '4')
    elif sort_len == 1:
        #last_type.append(motif_param)
        return motif_param
        exit()
    else:
        print("Error---check_motif")

    sort_all=sort_temp+sort_temp  
    sort_all_len=len(sort_all)
    count=sort_len
    
    for i in range(sort_len): 
        last_sort.append(sort_all[i:(sort_all_len - count)]) 
        count-=1
    sort_te = sorted(last_sort, key=lambda i: int(re.match(r'(\d+)', i).group()))  
    sort_te =sort_te[0] 
    sort_temp = "".join(sort_te) 

    if sort_len == 2 or sort_len == 4 or sort_len ==5 or sort_len ==3:  
        sort_temp=sort_temp.replace('1','A')
        sort_temp=sort_temp.replace('2','C')
        sort_temp=sort_temp.replace('3','G')
        sort_temp=sort_temp.replace('4','T')
    elif sort_len==6:
        sort_temp=sort_temp.replace('1','A')
        sort_temp=sort_temp.replace('2','C')
        sort_temp=sort_temp.replace('3','G')
        sort_temp=sort_temp.replace('4','T')
        sort_hexa = "TTAGGGTTAGGG"
        if sort_temp in sort_hexa:
            sort_temp ="TTAGGG"
    sort_temp_last = process_sequence_rev(sort_temp)
    return sort_temp_last
    
    
def txt_classify(chr_file):
    if '_all.txt' in chr_file:
        data_pd =pd.read_csv(in_outfile+chr_file,sep='\t',dtype={'    Iterations': np.str,'    Tract-size':np.str,'      Start': np.str,'      End': np.str})
        data_pd = data_pd.drop(0)
        data_pd = data_pd.drop(columns=[' Imperfection(p)%  '])
        data_pd = data_pd.reset_index(drop=True)
        data_pd['Consensus'] = data_pd['Consensus'].apply(lambda x: check_motif(x)) 
        data_pd.to_csv(in_outfile+chr_file.split('.txt')[0]+"-classify.txt",sep='\t',index =False)
        print(chr_file.split('.txt')[0]+"-classify.txt","txt type classification successful-6")
        return results_df

        

def txt_deal(num,filen):
    try:
        add=num
        name=filen
        with open(infile+filen+"_summary.txt","r",encoding="utf-8") as fp:
            data=fp.readlines()
        datas=data[2:]
        pf=open(outfile+filen+".txt","w",encoding="utf-8")
        pf.write(data[0])
        pf.write(data[1])
        for data in datas:
            data_item=re.findall(r'(\w+)(\s+)(\d+)(\s+)(\d+)(\s+)(\d+)(\s+)(\d+)(\s+)(\d+)',data)[0]
            first=str(int(data_item[6])+add)
            second=str(int(data_item[8])+add)
            pf.write(data_item[0]+data_item[1]+data_item[2]+data_item[3]+data_item[4]+data_item[5]+first+data_item[7]+second+data_item[9]+data_item[10]+"\n")
        pf.close()
        data_txt = pd.read_csv(outfile+filen+".txt", sep='\t')
        data_txt = data_txt.drop(0)
        data_txt = data_txt.reset_index(drop=True)
        data_list = data_txt.loc[data_txt['Consensus'].str.contains('N')] 
        # print()
        if len(data_list) != 0:
            print("-------------------------------------------------")
            print("\n"+infile+filen+"_summary.txt"+"------Data error, there is Gap data. Please check the original Fasta file")
            print(outfile+filen+".txt"+"------Data error, there is Gap data. Please check the original Fasta file"+"\n")
            print(data_list)
            print("\n")
            print("-------------------------------------------------")
        print("\n"+filen + "_summary.txt", " Txt processed successfully-1")
    except:
        print(filen+"_summary.txt","------Error")

# Core function

def ssr_density_statistics(filename, start, end, window, unit, segment):
    fileIn = codecs.open(fileInPath + filename, "r", encoding="utf-8")
    textIn = fileIn.readlines()[2:]
    fileIn.close()
    # unit_pre_window = int(window / unit)
    # fr_5w_num = (end - start) // window + 1
    new_start = start // unit * unit
    new_end = (end // unit + 1) * unit
    # fr_1_num = (new_end - new_start) // unit
    win_num = window // 1000
    fr_inx = list(range(new_start, end, unit))
    fr_inx.append(end)
    #print(fr_inx)
    # fr_inx = fr_inx[1:]
    inx = 0
    label = 0
    type_1 = segment
    data_ssr = []
    for i in range(1, len(fr_inx)):
        k_mono = k_di = k_tri = k_tetra = k_penta = k_hexa = k_mono_len = \
            k_di_len = k_tri_len = k_tetra_len = k_penta_len = k_hexa_len = 0
        pre_unit = fr_inx[i-1]+1 
        #print(pre_unit)
        curr_unit = fr_inx[i]
        #print(curr_unit)
        for k in textIn[inx:]:
            k_start = int(k.split()[3])
            k_end = int(k.split()[4])

            if k_start > pre_unit and k_end <= curr_unit:
                k_len = int(k.split()[2])
                k_type = len(k.split()[0])

                if k_type == 1:
                    k_mono += 1
                    k_mono_len += k_len
                    inx += 1
                elif k_type == 2:
                    k_di += 1
                    k_di_len += k_len
                    inx += 1
                    #print(k_di_len)
                elif k_type == 3:
                    k_tri += 1
                    k_tri_len += k_len
                    inx += 1
                elif k_type == 4:
                    k_tetra += 1
                    k_tetra_len += k_len
                    inx += 1
                elif k_type == 5:
                    k_penta += 1
                    k_penta_len += k_len
                    inx += 1
                elif k_type == 6:
                    k_hexa += 1
                    k_hexa_len += k_len
                    inx += 1
            elif k_start <= pre_unit and k_end <= curr_unit:
                ov_len = k_end - pre_unit + 1 
                k_type = len(k.split()[0])
                if k_type == 1:
                    k_mono += 1
                    k_mono_len += ov_len
                    inx += 1
                elif k_type == 2:
                    k_di += 1
                    k_di_len += ov_len
                    inx += 1
                elif k_type == 3:
                    k_tri += 1
                    k_tri_len += ov_len
                    inx += 1
                elif k_type == 4:
                    k_tetra += 1
                    k_tetra_len += ov_len
                    inx += 1
                elif k_type == 5:
                    k_penta += 1
                    k_penta_len += ov_len
                    inx += 1
                elif k_type == 6:
                    k_hexa += 1
                    k_hexa_len += ov_len
                    inx += 1
            elif k_end > curr_unit:
                if k_start <= curr_unit and pre_unit <= k_start:
                    #ov_len = k_end - curr_unit
                    k_len = curr_unit - k_start +1 
                    k_type = len(k.split()[0])
                    if k_type == 1:
                        k_mono += 1
                        k_mono_len += k_len
                    elif k_type == 2:
                        k_di += 1
                        k_di_len += k_len
                    elif k_type == 3:
                        k_tri += 1
                        k_tri_len += k_len
                    elif k_type == 4:
                        k_tetra += 1
                        k_tetra_len += k_len
                    elif k_type == 5:
                        k_penta += 1
                        k_penta_len += k_len
                    elif k_type == 6:
                        k_hexa += 1
                        k_hexa_len += k_len
                elif k_start <= curr_unit and k_start < pre_unit:
                    #ov_len = k_end - curr_unit
                    #len_temp = pre_unit - k_start
                    #ov_len = ov_len + len_temp
                    #k_len = int(k.split()[2]) - ov_len
                    k_len = curr_unit - pre_unit + 1 
                    k_type = len(k.split()[0])
                    if k_type == 1:
                        k_mono += 1
                        k_mono_len += k_len
                    elif k_type == 2:
                        k_di += 1
                        k_di_len += k_len
                    elif k_type == 3:
                        k_tri += 1
                        k_tri_len += k_len
                    elif k_type == 4:
                        k_tetra += 1
                        k_tetra_len += k_len
                    elif k_type == 5:
                        k_penta += 1
                        k_penta_len += k_len
                    elif k_type == 6:
                        k_hexa += 1
                        k_hexa_len += k_len
                break
            #if 
        k_total = k_mono + k_di + k_tri + k_tetra + k_penta + k_hexa
        k_total_len = k_mono_len + k_di_len + k_tri_len + k_tetra_len + \
                      k_penta_len + k_hexa_len

        k_position = fr_inx[i]
        k_label = segment + "-Z" + str(label//win_num + 1)
        label += 1
        ssr_pre = [type_1, k_position, k_label, k_mono, k_di, k_tri, k_tetra, k_penta, k_hexa, k_total,
                  k_mono_len, k_di_len, k_tri_len, k_tetra_len, k_penta_len, k_hexa_len, k_total_len,
                  ex_num, ex_den]
        if k_total_len >1000:
              print("total_len Total length greater than 1000"+"------"+str(k_position)+"------"+str(k_label)+"\n")
        data_ssr.append(ssr_pre)

    fileOutname = ssr_allPath + filename.split(".")[0] + "_all.csv"
    # print(fileOutname)
    textOut = codecs.open(fileOutname, 'w')
    textOut.write(fileHead)
    for ssr in data_ssr:
        fi = [str(x) for x in ssr]
        textOut.write(','.join(fi) + '\n')
    textOut.close()
    print("\n"+filename.split(".")[0] + "_all.csv"+ "  CSV generation-2")
    

def csv_deal(name_c,map_add,flag_head):
    name_all = name_c.split("-S")[0] + "_all.csv"
    name_c = all_file +name_c
    data = pd.read_csv(name_c)
    count_c = data['label'].value_counts(ascending=True)
    #print(count_c[0])#print(count_c.index[0])#print(len(count_c))#print(count_c.index)
    if count_c[0] < flag_num:  
        map_num = count_c.index[0].split('-Z')
        map_re = str(int(map_num[1])-1)
        if int(map_num[1]) != 1:
          new_map = map_num[0]+'-Z'+map_re
            #print(new_map)
          data = data.replace(count_c.index[0],new_map)
    data['label_2'] = data["label"].map(lambda x: int(x.split('-Z')[1]))
    #print(data['label_2'])
    data['label_2'] = data['label_2'] +map_add
    #print(data['label_2'])
    data["label_2"]= data["label_2"].map(lambda x: str(x))
    data['label'] = data["label"].map(lambda x: x.split('-Z')[0])
    data['label'] = data['label']+'-Z'+data['label_2']
    #print(data['label'])
    del data['label_2']
    #print(data['label'].iloc[-1].split('-Z')[1])
    now_map = int(data['label'].iloc[-1].split('-Z')[1])
    if flag_head == 0:
        data.to_csv(all_file+ALL_mer+name_all,sep=',',index = False,mode ='a+')
    else:
        data.to_csv(all_file+ALL_mer+name_all,sep=',',index = False,mode ='a+',header = False)
    #print(data)
    #print(now_map)
    return now_map

def txt_merge(txt_name,txt_cm):
    with open(outfile+txt_name+".txt","r",encoding="utf-8") as fp:
        #fp.readline()
        data=fp.readlines()
    datas=data[2:]
    with open('in_outfile\\'+txt_name.split('-S')[0]+"_all.txt","a+",encoding="utf-8") as pf:
        if txt_cm == 0:
            pf.write(data[0])
            pf.write(data[1])
        for temp_data in datas:
            pf.write(temp_data)
        #pf.write(datas)
    #print('%d: '%txt_cm,txt_name+".txt","Successfully processed")    
    

    
    
def Landscape_bedgf(i_file):
    if len(seg_list) != 0:
        global data_seg
        gap_flag = '1'
    else:
        gap_flag = '0'
    data_pd = pd.read_csv(mer_file+i_file)
    data_pd=data_pd[['position','ratio_len_total']]
    data_pd = data_pd.reset_index(drop=True)  
    chr_temp = i_file.split('_all')[0]
    if gap_flag == '1':
        data_seg_templ = data_seg[data_seg['filename'].str.contains(chr_temp)]
    data_pd['s_position'] = data_pd['position'] - 1000  
    data_pd = data_pd.astype({'s_position': 'str'})
    data_pd['check_end'] = data_pd['s_position'].str.endswith('000')  
    data_pd.loc[(data_pd['s_position'] == '0'),'check_end'] = 'True'  
    list_se = data_pd[(data_pd['check_end'] == False)].index.tolist() 
    for i_e in list_se: 
        s_temp = data_pd.loc[i_e]['s_position']
        npeak_e = str(s_temp)[-3:]
        e_len = len(str(s_temp))
        lnpeak_e = str(s_temp)[:-3]
        s_position = int(lnpeak_e.ljust(e_len, '0'))  
        s_position = str(s_position+1000)
        data_pd.loc[i_e, ('s_position')] = s_position 
    data_pd = data_pd.astype({'s_position': 'int'})
    #print(data_pd)
    if gap_flag == '1':
        for rows in data_seg_templ.index: 
            seg_data = data_seg_templ.loc[rows]['from']
            data_pd['Diff_value'] = seg_data - data_pd['s_position']  
            data_pd['check_value']  = (data_pd["Diff_value"]>1)&(data_pd['Diff_value']<1001)  
            list_sv = data_pd[(data_pd['check_value'] == True)].index.tolist() 
            #print(list_sv)
            if len(list_sv) != 0: 
                for i_sv in list_sv:
                    #print(data_pd)
                    if data_pd.loc[i_sv, ('position')] > seg_data:
                        data_pd.loc[i_sv, ('s_position')] = seg_data -1  
                        #print(seg_data -1)
            data_pd = data_pd.drop(columns=['Diff_value','check_value'])
    #print(data_pd)
    peak_n = chr_temp.split('-')
    peak_chr = peak_n[0]+'-'+peak_n[1]
    if '0' in peak_n[1] and peak_n[1] != '10' and peak_n[1] != '20':
        chr_name = 'chr'+peak_n[1].split('0')[1]
    else:
        chr_name = 'chr' + peak_n[1]
    data_pd['chr_name'] =chr_name
    order = ['chr_name','s_position','position','ratio_len_total']
    data_pd = data_pd[order]
    #print(data_pd)
    data_pd.to_excel(land_bgf+i_file.split('.csv')[0]+'_bedgraph.xlsx', index=False,engine="openpyxl")
    print(i_file.split('.csv')[0]+'_bedgraph.xlsx'+"------"+"generated")




txt_dc = config_control[0]["txt_deal"]   #Txt processing
csv_gc = config_control[0]["csv_gener"]  #Generate CSV
csv_mc = config_control[0]["csv_mer"]    #CSV Merge
txt_mc = config_control[0]["txt_mer"]   #Merge txt
txt_cc = config_control[0]["txt_classify"]   #Txt_ Motif classification
land_bc = config_control[2]["land_bedg"] #Used to control whether to generate a bedgraph file



num_cores = multiprocessing.cpu_count()-4

if txt_dc == '1':
    Parallel(n_jobs=num_cores)(delayed(txt_deal)(temp["num"],temp["filename"]) for temp in config_txt)
    print("Txt processing completed-1")

print("-------------------------------------------------------"+"\n")

if csv_gc == '1':
    Parallel(n_jobs=num_cores)(delayed(ssr_density_statistics)(item["name"]+".txt",item["start"],item["end"],50000,1000,item["segname"]) for item in config)
    print("CSV generation-2")

print("-------------------------------------------------------"+"\n")


#name_temp = config[0]['name'].split('-S')[0]

if csv_mc == '1':
    for item in config_peak:
        name_temp = item['name']
        file_path = all_file+ALL_mer+name_temp+"_all.csv"
        if (os.path.exists(file_path)):
            os.remove(file_path)
if csv_mc == '1':
    for item in config:   
        csv_name = item["name"]+"_all.csv"
        csv_chr = item["name"].split('-S')[0]
        if count_csv == 0:
            old_csv = csv_chr
        if csv_chr != old_csv:
            count_csv = 0
            old_csv = csv_chr
            map_add =0
        map_add = csv_deal(csv_name,map_add,count_csv)
        print("%s"%csv_name+"  "+"CSV fragments merged successfully-3")
        count_csv += 1
print("-------------------------------------------------------"+"\n")



if txt_mc == '1':
    for item in config_peak:
        name_temp = item['name']
        file_txt = in_outfile+name_temp+"_all.txt"
        #print(file_txt)
        if (os.path.exists(file_txt)):
            os.remove(file_txt)

if txt_mc == '1': 
    for temp_i in config:
        txt_n = temp_i['name']
        txt_chr = txt_n.split('-S')[0]
        if txt_c == 0:
            old_chr = txt_chr
        if txt_chr != old_chr:
            txt_c = 0
            old_chr = txt_chr
        txt_merge(txt_n,txt_c)
        txt_c += 1
        print('%d: '%txt_c,txt_n+".txt","Txt merged successfully-4")


print("-------------------------------------------------------"+"\n")

#(columns=['    Iterations','      Start','      End',' Imperfection(p)%  '])
if txt_cc == '1':
    chr_files = os.listdir(in_outfile)
    #print(chr_files)
    results_data = Parallel(n_jobs=num_cores)(delayed(txt_classify)(chr_file) for chr_file in chr_files)
    results_dfs = pd.concat(results_data)
    print("\n"+"Txt type classification successful-5")
    print(results_dfs)

    results_dfsl = results_dfs.drop_duplicates(subset='Orig_Seq', keep='first')
    results_dfslc = results_dfsl.copy()
    results_dfslc['len'] = results_dfslc['Orig_Seq'].str.len()
    results_dfslsl = results_dfslc.sort_values(by=['len', 'Orig_Seq'])

    results_dfslsls = results_dfslsl.drop_duplicates(subset='Orig_Seq', keep='first')
    results_dfslsls.to_excel('output.xlsx', index=False)
print("-------------------------------------------------------"+"\n")



count_bg = 0 #Used to merge bedG files
if land_bc == '1':
    deal_csvl = os.listdir(mer_file)  #Get mer_file name under file is returned as a list
    Parallel(n_jobs=num_cores)(delayed(Landscape_bedgf)(csv_f) for csv_f in deal_csvl)
    bedGf = os.listdir(land_bgf)      #Used for loop merging files
    name_bedG = bedGf[0].split('-')[0] +"-landscape_all"+".bedGraph"
    for i_f in bedGf:
        data_bedG = pd.read_excel(land_bgf + i_f,engine="openpyxl")
        if count_bg == 0:
            if (os.path.exists(bed_outfile+name_bedG)):
                os.remove(bed_outfile + name_bedG)
            data_bedG.to_csv(bed_outfile + name_bedG, sep='\t', index=False, mode='a+')
        else:
            data_bedG.to_csv(bed_outfile+ name_bedG, sep='\t', index=False, mode='a+', header=False)
        count_bg += 1
    print("\n"+"xDMA_ucsc BedGraph file completed---6")
    print("-----------------------------------------------------------------------------------------------------")
    print("\n")


print("-------------------------------------------------------"+"\n")
print("Processing completed ----------------------Software version：V3.0")
input()