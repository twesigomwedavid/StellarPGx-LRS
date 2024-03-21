#!/usr/bin/env python3

import os
import sys
import math
        

def get_total_CN(cov_file):

    all_reg = []
    for line in open(cov_file, "r"):
        line = line.strip().split()
        all_reg.append(line)

    av_2d6_cov = float(all_reg[2][3])/(float(all_reg[2][2]) - float(all_reg[2][1]))
    av_in1_3pr = float(all_reg[1][3])/(float(all_reg[1][2]) - float(all_reg[1][1]))
    av_ex9_3pr = float(all_reg[0][3])/(float(all_reg[0][2]) - float(all_reg[0][1]))
    av_ctrl_reg = float(all_reg[3][3])/(float(all_reg[3][2]) - float(all_reg[3][1]))

    # av_ctrl_cov = (av_3a5_cov + av_3a4_cov)/2
    av_ctrl_cov = av_ctrl_reg
    
    comp_av = av_2d6_cov/av_ctrl_cov
    temp_cn = 2 * comp_av
    total_cn = round(temp_cn)
    
    in1_3pr = round(2 * av_in1_3pr/av_ctrl_cov) 
    ex9_3pr = (2 * av_ex9_3pr/av_ctrl_cov)

    # return [round(av_2d6_cov), (av_ctrl_cov)]
    
    return [str(int(total_cn)), round(av_2d6_cov), str(int(in1_3pr)), round(av_ctrl_cov), str(ex9_3pr), round(av_in1_3pr)];


samp_gt = ""
samp_gt_hap1 = ""

# def del_test(sv_del):

#     if os.stat(sv_del).st_size == 0:
#         return "None"

#     else:
#         for line in open(sv_del, "r"):
#             if "COVERAGE" in line:
#                 line = line.strip().split()
#                 ABHom = line[-1]
#                 ABHet = line[-2]
#                 GT = line[2]
#                 DP = int(line[3])
        
#                 if float(ABHom) == 1.0:
#                     return "*5/*5"
#                 elif float(ABHom) == -1.0:
#                     return "*5"
#             else:
#                 pass


# hap_adv_list = []
# hap_t1 = []


# def del_adv_test(hap_dbs, cand_allele1, cand_allele2, test_allele1, test_allele2, core_vars):
#     g = open(hap_dbs, "r")
#     for line in g:
#         line = line.strip().split()
#         hap_adv_list.append(line)
    
#     a1 = core_vars.split(";")

#     for i in a1:
#         if i[-3:] == "0/1":
#             hap_t1.append(i[:-4])

  
#     for elem in hap_adv_list:
#         if elem[1] == cand_allele1: 
#             list_t1 = (elem[2]).split(';')

#         if elem[1] == cand_allele2:
#             list_t2 = (elem[2]).split(';')
    
#     if hap_t1[0] in list_t1:
#         return test_allele1

#     elif hap_t1[0] in list_t2:
#         return test_allele2


het_hom_list = []
het_hom_list_new = []

def dup_test_init(sv_dup, av_cov):
    for line in open(sv_dup, "r"):
    #     if "COVERAGE" in line:
    #         continue
    #     elif "AGGREGATED" in line:
    #         continue

    #     else:
        fields = line.strip().split()
        het_hom_list.append(fields)

    test_list1 = []

    for i in het_hom_list:
        test_list1.append(int(i[2]))

    av_read_cov = sum(test_list1)/len(test_list1)
    norm_cov = (av_cov + av_read_cov)/2

    for i in het_hom_list:
        supp_reads = round(float(i[-2])*int(i[2]))
        i.append(round(supp_reads/norm_cov, 3))
        i.append(supp_reads)
        het_hom_list_new.append(i)


    return (het_hom_list_new)


dp_list = []

def test_base(sv_dup, av_cov, supp_core_vars1, supp_core_vars2, cand_allele1, cand_allele2):
    
    for line in open(sv_dup, "r"):
        line = line.strip().split()
        dp_list.append(line)

    het_list1 = []

    het_list2 = []

    het_list3 = []

    abhet_list = []

    abhet_list2 = []

    # ad_ratio_list = []

    # ad_ratio_list2 = []

    
    if supp_core_vars1 == "ref_haplotype":
        supp_core_vars1 = []
    else:
        supp_core_vars1 = supp_core_vars1.split(";")

    if supp_core_vars2 == "ref_haplotype":
        supp_core_vars2 = []
    else:
        supp_core_vars2 = supp_core_vars2.split(";")
    
    all_core_vars = supp_core_vars1 + supp_core_vars2 

    
    for i in dp_list:
        
        if (i[2] == "0|1" or i[2] == "1|0") and i[-1] != ".":
            het_list1.append(i[0])
            ad_test = i[4].split(",") 
            abhet = round((int(ad_test[1])/int(ad_test[0]) + int(ad_test[1])), 4)

            abhet_list.append(abhet)

        
        if i[0] in all_core_vars and (i[2] == "0|1" or i[2] == "1|0") and i[-1] != ".":
            ind1 = het_list1.index(i[0])
            abhet_list2.append(abhet_list[ind1])

        else:
            pass
            
    max_het = max(abhet_list2)
    max_het_pos = abhet_list2.index(max_het)
    var = all_core_vars[max_het_pos]
    # var2 = all_core_vars[max_ratio_pos]

    return [var, max_het, supp_core_vars1, supp_core_vars2];


def test_dup1(sv_dup, av_cov, supp_core_vars1, supp_core_vars2, cand_allele1, cand_allele2, cn):
    c_num = ""
    in_list = test_base(sv_dup, av_cov, supp_core_vars1, supp_core_vars2, cand_allele1, cand_allele2)

    var = in_list[0]
    max_het = in_list[1]
    # var2 = in_list[2]
    # max_ratio = in_list[3]
    supp_core_vars1 = in_list[2]
    supp_core_vars2 = in_list[3]


    res_dip = cand_allele1 + "/" + cand_allele2

    c_num = 2
    
    # if round(max_ratio) == 2:
    #     if var2 in supp_core_vars1:
    #         res_dip = cand_allele2 + "/" + cand_allele1 + "x" + "2"
    #         c_num = 3
            
    #     elif var2 in supp_core_vars2:
    #         res_dip = cand_allele1 + "/" + cand_allele2 + "x" + "2"
    #         c_num = 3
            
    # else:
    #     res_dip = cand_allele1 + "/" + cand_allele2
    #     c_num = 2
        
    return [res_dip, c_num];    
            

def test_dup2(sv_dup, av_cov, supp_core_vars1, supp_core_vars2, cand_allele1, cand_allele2, cn):
    c_num = ""
    in_list = test_base(sv_dup, av_cov, supp_core_vars1, supp_core_vars2, cand_allele1, cand_allele2)

    var = in_list[0]
    max_het = in_list[1]
    max_hap_cn = round(max_het * int(cn))
    hap2_cn = int(cn) - max_hap_cn 
    supp_core_vars1 = in_list[2]
    supp_core_vars2 = in_list[3]


    if var in supp_core_vars1:
        if hap2_cn == 1:
            res_dip = cand_allele2 + "/" + cand_allele1 + "x" + max_hap_cn

        elif hap2_cn != 1:
            res_dip = cand_allele2 + "x" + hap2_cn + "/" + cand_allele1 + "x" +	max_hap_cn

    elif var in supp_core_vars2:
        if hap2_cn == 1:
            res_dip = cand_allele1 + "/" + cand_allele2 + "x" + max_hap_cn

        elif hap2_cn !=	1:
            res_dip = cand_allele1 + "x" + hap2_cn + "/" + cand_allele2 + "x" +	max_hap_cn
            
    return res_dip



hap_def_list = []
allele_cn_list = []

def dup_test_cn_3_4(sv_dup, hap_dbs, cand_allele1, cand_allele2, test_allele1, test_allele2, c_num, av_cov, in_list):

    g = open(hap_dbs, "r")
    for line in g:
        line = line.strip().split()
        hap_def_list.append(line)
        
            
    test_list1 = []
    test_list2 = []
    het_list = []


    for i in in_list:
        if i[1] == "0/1":
            het_list.append(i)

    for i in het_list:
        test_list1.append(i[0])
        test_list2.append(i[-2])


    max_het = max(test_list2)
    max_het_pos = test_list2.index(max_het)
    var = test_list1[max_het_pos]
    
    for elem in hap_def_list:
        if elem[1] == cand_allele1:
            list_3t = elem
            list_3t_2 = list_3t[2].split(';')
            l3 = len(list_3t_2)
            
        if elem[1] == cand_allele2:
            list_4t = elem
            list_4t_2 = list_4t[2].split(';')
            l4 = len(list_4t_2)

    hdb_list = list_3t_2 + list_4t_2


    index_var = hdb_list.index(var)

    if index_var < l3:
        allele_cn_list.append(test_allele1)
        allele_cn_list.append(int(round(max_het*int(c_num))))

    elif index_var >= l3:
        allele_cn_list.append(test_allele2)
        allele_cn_list.append(int(round(max_het*int(c_num))))

    
    if allele_cn_list[0] == test_allele1:
        rt_2 = int(c_num) - allele_cn_list[1]
        allele_cn_list.append(test_allele2)
        allele_cn_list.append(rt_2)

    elif allele_cn_list[0] == test_allele2:
        rt_2 = int(c_num) - allele_cn_list[1]
        allele_cn_list.append(test_allele1)
        allele_cn_list.append(rt_2)


    if allele_cn_list[1] == 0:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3] - 1)

    elif allele_cn_list[3] == 0:
        res_dip = allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1] - 1)

    elif allele_cn_list[1] == 1:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3])

    elif allele_cn_list[3] == 1:
        res_dip = allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1])  

    elif allele_cn_list[1] == 2:
        res_dip = allele_cn_list[0] + "x2" + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3])

    elif allele_cn_list[3] == 2:
        res_dip = allele_cn_list[2] + "x2" + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1])

    elif allele_cn_list[1] == -1:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3] - 2)

    elif allele_cn_list[3] == -1:
        res_dip =  allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1] - 2)

    elif allele_cn_list[3] < -1:
        res_dip =  allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(3)

    elif allele_cn_list[1] < -1:
        res_dip =  allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(3)


    else:            
        res_dip = 'check'

    return res_dip




def dup_test_cn_n(sv_dup, hap_dbs, cand_allele1, cand_allele2, test_allele1, test_allele2, c_num, av_cov, in_list):

    g = open(hap_dbs, "r")
    for line in g:
        line = line.strip().split()
        hap_def_list.append(line)


    test_list1 = []
    test_list2 = []
    het_list = []


    for i in in_list:
        if i[1] == "0/1":
            het_list.append(i)

    for i in het_list:
        test_list1.append(i[0])
        test_list2.append(i[-2])

    max_het = max(test_list2)
    max_het_pos = test_list2.index(max_het)
    var = test_list1[max_het_pos]



    for elem in hap_def_list:
        if elem[1] == cand_allele1:
            list_3t = elem
            list_3t_2 = list_3t[2].split(';')
            l3 = len(list_3t_2)

        if elem[1] == cand_allele2:
            list_4t = elem
            list_4t_2 = list_4t[2].split(';')
            l4 = len(list_4t_2)

    hdb_list = list_3t_2 + list_4t_2


    index_var = hdb_list.index(var)

    if index_var < l3:
        allele_cn_list.append(test_allele1)
        allele_cn_list.append(int(round(max_het*int(c_num)-0.15)))

    elif index_var >= l3:
        allele_cn_list.append(test_allele2)
        allele_cn_list.append(int(round(max_het*int(c_num)-0.15)))


    if allele_cn_list[0] == test_allele1:
        rt_2 = int(c_num) - allele_cn_list[1]
        allele_cn_list.append(test_allele2)
        allele_cn_list.append(rt_2)

    elif allele_cn_list[0] == test_allele2:
        rt_2 = int(c_num) - allele_cn_list[1]
        allele_cn_list.append(test_allele1)
        allele_cn_list.append(rt_2)

    if allele_cn_list[1] == 0:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3] - 1)

    elif allele_cn_list[3] == 0:
        res_dip = allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1] - 1)

    elif allele_cn_list[1] == 1:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3])

    elif allele_cn_list[3] == 1:
        res_dip = allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1])

    elif allele_cn_list[1] == 2:
        res_dip = allele_cn_list[0] + "x2" + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3])

    elif allele_cn_list[3] == 2:
        res_dip = allele_cn_list[2] + "x2" + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1])

    elif allele_cn_list[1] == 3:
        res_dip = allele_cn_list[0] + "x3" + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3])

    elif allele_cn_list[3] == 3:
        res_dip = allele_cn_list[2] + "x3" + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1])

    elif allele_cn_list[1] == 4:
        res_dip = allele_cn_list[0] + "x4" + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3])

    elif allele_cn_list[3] == 4:
        res_dip = allele_cn_list[2] + "x4" + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1])

    elif allele_cn_list[1] == -1:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3] - 2)

    elif allele_cn_list[3] == -1:
        res_dip =  allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1] - 2)

    elif allele_cn_list[3] < 0:
        res_dip =  allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1] + allele_cn_list[3] - 1)

    elif allele_cn_list[1] < 0:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3] + allele_cn_list[1] - 1)
        

    else:
       res_dip = 'check'

    return res_dip



def hybrid_test_68(sv_dup, c_num, av_cov, cn_in1_3pr1, in_list):
    
    test_list1 = []
    test_list2 = []
    test_list3 = []


    for i in in_list:
        test_list1.append(i[0])
        test_list2.append(abs(float(i[-2])))
        test_list3.append(i[-1])

    index1 = test_list1.index('42130692~G>A')
    index2 = test_list1.index('42128945~C>T')
    
    val_68 = test_list3[index1]
    val_4 = test_list3[index2]
    
    rt = val_68/val_4

    
    if rt <= 1.4:
        return 'norm_dup'

    elif rt > 1.4:
        return 'hyb_68'

       
    else:
        return 'norm_dup'


def hyb_test_5_68_4(sv_del, in1_3pr1_float, av_cov):
    test_del = []
    for line in open(sv_del, "r"):
        if "COVERAGE" in line:
            test_del.append(line.strip())

    # if len(test_del) == 0:
    #     return 'norm_art'

    # elif len(test_del) > 0:
    #     return 'del_hyb'

    t1 = 2 * in1_3pr1_float/av_cov

    if len(test_del) == 0 and (1.6 < t1 < 2.8):
        return 'norm_art'

    elif len(test_del) > 0 and t1 < 1.6:
        return 'del_hyb'


def hybrid_test_36(sv_dup, cn, av_cov, cn_ex9_3pr, cn_2d7_ex9, cn_2d7_in4_in8):


    if ((int(cn) - 1) - 0.3) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5) or (2.5 < (2 * float(cn_2d7_ex9)/float(cn_2d7_in4_in8)) < 3.5):
        return 'hyb_36_10'

    elif (int(cn) - 2) <= float(cn_ex9_3pr) < (int(cn) - 2 + 0.7):
        return 'hyb_36_36'
    else:
        return 'norm_dup'


def hybrid_test_36_single(sv_dup, cn, av_cov, cn_ex9_3pr):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_star10'

    elif ((int(cn) - 1) - 0.35) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        return 'hyb_36_single'

    else:
        return 'norm_star10'



def hybrid_test_36_mod(sv_dup, cn, av_cov, cn_ex9_3pr):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_mt'
                                                                                                  
    elif ((int(cn) - 1) - 0.05) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        return 'hyb_36_10'
    elif (int(cn) - 2) <= float(cn_ex9_3pr) < (int(cn) - 2 + 0.95):
        return 'hyb_36_36'


def hybrid_test_36_multi(sv_dup, cn, av_cov, cn_ex9_3pr):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_mt'

    elif ((int(cn) - 1) - 0.05) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        return 'hyb_36_10'
    elif (int(cn) - 2) <= float(cn_ex9_3pr) < (int(cn) - 2 + 0.95):
        return 'hyb_36_36'
    elif (int(cn) - 3) <= float(cn_ex9_3pr) < (int(cn) - 3 + 0.95):
        return 'hyb_36_36_36'    
    else:
        return 'check'


def hybrid_test_36_multi_10(sv_dup, cn, av_cov, cn_ex9_3pr, cn_star10):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_mt'

    elif float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        cn_star36 = int(cn) - int(round(float(cn_ex9_3pr)))
        adj_cn_star10 = int(cn_star10) - cn_star36
        
        if cn_star36 == 1:
            return '*36+*10x' + str(adj_cn_star10)
        else:
            return '*36x' + str(cn_star36) + '+*10x'  + str(adj_cn_star10) 

    else:
        return 'check'


def hybrid_13_2_v1(cov_in4_3pr, cov_5pr_in4):
    
    if 0.85 < float(cov_in4_3pr)/float(cov_5pr_in4) < 1.2:
        return 'norm_var'
    elif 0.45 < float(cov_in4_3pr)/float(cov_5pr_in4) < 0.75:
        return 'hyb_13_2'
    else:
        return 'norm_var'


def hybrid_13_2_v2(cov_2d7_ex2_in8, cov_2d7_5pr_in1):

    if 0.85 < float(cov_2d7_ex2_in8)/float(cov_2d7_5pr_in1) < 1.2:
        return 'norm_var'
    elif 0.45 < float(cov_2d7_ex2_in8)/float(cov_2d7_5pr_in1) < 0.75:
        return 'hyb_13_2_v2'
    else:
        return 'norm_var'



def tandem_90_1(in_list, alt_allele, cn):

    test_list1 = []
    test_list2 = []
    test_list3 = []
                                                                                                                                     
    for i in in_list:
        test_list1.append(i[0])
        test_list2.append(abs(float(i[-2])))
        test_list3.append(i[-1])


    if len(test_list1) > 1:
        index1 = test_list1.index('42129098~T>C')
        a = test_list3[index1]
        test_list3.pop(index1)
        b = max(test_list3)

        c = round(b/a)
        
        if int(cn) == 3 and c == 1:
            res = alt_allele + "/" + "*90+*1"

        elif int(cn) == 3 and c > 1:
            res = alt_allele + "x2" + "/" + "*90"

        elif int(cn) == 4 and c == 2:
            res = alt_allele + "x2" + "/" + "*90+*1"
    
        elif int(cn) == 4 and c >= 3:
            res = alt_allele + "x3" + "/" + "*90"

    else:
        val1 = test_list2[0]
        val2 = round(val1 * int(cn))

        if int(cn) == 3 and val2 == 1:
            res = '*1/*90+*1'
        elif int(cn) == 3 and val2 == 2:
            res = '*90/*90+*1'
        elif int(cn) == 4 and val2 == 1:
            res = '*1x2/*90+*1'
        elif int(cn) == 4 and val2 == 2:
            res = '*90+*1/*90+*1'
        elif int(cn) == 4 and val2 == 3:
            res = '*90x2/*90+*1'
        elif int(cn) == 3 and val2 == 3:
            res = '*90/*90x2'
        elif int(cn) == 4 and val2 == 4:
            res = '*90/*90x3'

    return res


def tandem_57_10(in_list, alt_allele, cn):

    test_list1 = []
    test_list2 = []
    test_list3 = []

    for i in in_list:
        test_list1.append(i[0])
        test_list2.append(abs(float(i[-2])))
        test_list3.append(i[-1])

    if len(test_list1) > 1:
        index1 = test_list1.index('42129906~G>A')
        a = test_list3[index1]
        test_list3.pop(index1)
        index2 = test_list1.index('42130692~G>A')
        m = test_list3[index2]
        test_list3.pop(index2)
        b = max(test_list3)

        c = round(b/a)
        p = round(m/a)

        if int(cn) == 3 and c == 1 and p > 1:
            res = alt_allele + "/" + "*57+*10"

        elif int(cn) == 3 and c > 1 and p == 1:
            res = alt_allele + "x2" + "/" + "*57"

        elif int(cn) == 4 and c == 2 and p > 1:
            res = alt_allele + "x2" + "/" + "*57+*10"

        elif int(cn) == 4 and c >= 3 and p == 1:
            res = alt_allele + "x3" + "/" + "*57"

        elif int(cn) == 4 and p == 1 and alt_allele == '*10':
            res = "*57+*10" + "/" + "*57+*10"

        elif int(cn) == 4 and p > 1 and alt_allele == '*10':
            res = "*10x2" + "/" + "*57+*10"

        else:
            res = alt_allele + "/" + "*57+*10"

    return res



def hybrid_test_83_single(sv_dup, cn, av_cov, cn_ex9_3pr):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_star39'

    elif ((int(cn) - 1) - 0.35) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        return 'hyb_83_single'

    else:
        return 'norm_star39'


def hybrid_test_83(sv_dup, cn, av_cov, cn_ex9_3pr):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_star39'

    elif ((int(cn) - 1) - 0.35) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        return 'hyb_83'

    else:
        return 'norm_star39'
