# bioinformatics

Команда
- [Лебедев Аркадий ИУ9-71Б](https://github.com/ark2016)
- [Кочетков Денис ИУ9-71Б](https://github.com/HumsterProgrammer)

## Описание задачи

Дано: файл с набором белков и образец белка

Найти: топ N похожих на образец белков из файла

Степень похожести должна определяться с учетом количества измений. Изменения бывают атомарные, то есть добавляющие, изменяющие или удаляющие один символ, или групповые. Групповые изменения также включают перетасовку групп: например, были группы ABC стали ACB. Чем больше изменений, тем менее похожими должны быть белки.

## Алгоритм реализации

Алгоритм основан на идее сжатия строки с помощью клеточного автомата (КА) и последующего сравнения строк с помощью оценки косинуса угла между векторными представлениями строк.

Алгоритм работы можно представить следующим образом:
1. Парсится файл и выделяются строки и подписи к ним.
2. Строки преобразуются к внутреннему представлению КА.
3. Далее, для всех строк из файла производится:
	- приведение строк к одной длинне с помощью КА
	- вычисляется косинус угла с помощью скалярных произведений
4. Производится сортировка по косинусу и выводится топ похожих белков. По-умолчанию, топ 100.

### Оптимизации

Алгоритм оптимизирован с помощью векторизации вычислений, паралельных вычислений, кэширования, препроцессинга etc.

## Тестирование
Было проведено тестирование на наборах данных. Здесь пока пусто.

### uniprot_sprot.fasta
```
(venv) vudrav@vudrav-XL442:~/Документы/0_sem7/bio/lab1/bioinformatics$ python3 read.py
Loaded 573661 records
                     id                                        description                                           sequence
0  sp|Q6GZX4|001R_FRG3G  sp|Q6GZX4|001R_FRG3G Putative transcription fa...  MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQV...
1  sp|Q6GZX3|002L_FRG3G  sp|Q6GZX3|002L_FRG3G Uncharacterized protein 0...  MSIIGATRLQNDKSDTYSAGPCYAGGCSAFTPRGTCGKDWDLGEQT...
2   sp|Q197F8|002R_IIV3  sp|Q197F8|002R_IIV3 Uncharacterized protein 00...  MASNTVSAQGGSNRPVRDFSNIQDVAQFLLFDPIWNEQPGSIVPWK...
3   sp|Q197F7|003L_IIV3  sp|Q197F7|003L_IIV3 Uncharacterized protein 00...  MYQAINPCPQSWYGSPQLEREIVCKMSGAPHYPNYYPVHPNALGGA...
4  sp|Q6GZX2|003R_FRG3G  sp|Q6GZX2|003R_FRG3G Uncharacterized protein 3...  MARPLLGKTSSVRRRLESLSACSIFFFLRKFCQKMASLVFLNSPVY...
Data saved as .\data\uniprot_data.csv
```

```
(venv) vudrav@vudrav-XL442:~/Документы/0_sem7/bio/lab1/bioinformatics$ python main_numba.py "ASDIFSDUYFGSUDYFGSUDYGFUSDY"
Searching for: ASDIFSDUYFGSUDYFGSUDYGFUSDY
In file: data/uniprot_data.csv
Используется CPU оптимизация (12 потоков)
Loading and parsing genome data...
Comparing genomes...
Comparing genomes: 100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 573661/573661 [00:43<00:00, 13294.52it/s]
Finding similar genomes...
Filtering results: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 573661/573661 [00:00<00:00, 2893609.48it/s]

Found 100 similar genomes:
sp|P05487|CONO_CONST Conopressin-S OS=Conus striatus OX=6493 PE=1 SV=1: 0.9128
sp|P16339|DNF1_LOCMI Locupressin OS=Locusta migratoria OX=7004 PE=1 SV=1: 0.8975
sp|P85308|LAC1B_CERUI Laccase-1b (Fragment) OS=Cerrena unicolor OX=90312 PE=1 SV=1: 0.8445
sp|P83279|FAR6_MACRS FMRFamide-like neuropeptide FLP6 OS=Macrobrachium rosenbergii OX=79674 PE=1 SV=1: 0.8425
sp|Q96DX7|TRI44_HUMAN Tripartite motif-containing protein 44 OS=Homo sapiens OX=9606 GN=TRIM44 PE=1 SV=1: 0.8272
sp|P83570|GWA_SEPOF Neuropeptide GWa OS=Sepia officinalis OX=6610 PE=1 SV=1: 0.8109
sp|P82099|EI03_LITRU Electrin-3 OS=Litoria rubella OX=104895 PE=1 SV=1: 0.8056
sp|P0DV81|BPPL_BITRH Bradykinin-potentiating peptide-like OS=Bitis rhinoceros OX=715877 PE=1 SV=1: 0.8054
sp|B3A0A2|FAR3_AUSRA Extended FMRFamide-3 OS=Austrophasma rawsonvillense OX=253137 PE=1 SV=1: 0.7957
sp|B3A0E1|FAR3_AUSGA Extended FMRFamide-3 OS=Austrophasma gansbaaiense OX=253136 PE=1 SV=1: 0.7871
sp|B3A064|FAR3_KARBI Extended FMRFamide-3 OS=Karoophasma biedouwense OX=253133 PE=1 SV=1: 0.7871
sp|B0M8U2|FAR3_KARBO Extended FMRFamide-3 OS=Karoophasma botterkloofense OX=253132 PE=1 SV=1: 0.7871
sp|P86088|AMY2_CAPCH Alpha-amylase 2 (Fragment) OS=Capsicum chinense OX=80379 PE=1 SV=1: 0.7856
sp|P58707|FLRN_ANTEL Antho-RNamide OS=Anthopleura elegantissima OX=6110 PE=1 SV=1: 0.7841
sp|Q03E51|RL7_PEDPA Large ribosomal subunit protein bL12 OS=Pediococcus pentosaceus (strain ATCC 25745 / CCUG 21536 / LMG 10740 / 183-1w) OX=278197 GN=rplL PE=3 SV=1: 0.7676
sp|A5F2X6|MIAB_VIBC3 tRNA-2-methylthio-N(6)-dimethylallyladenosine synthase OS=Vibrio cholerae serotype O1 (strain ATCC 39541 / Classical Ogawa 395 / O395) OX=345073 GN=miaB PE=3 SV=1: 0.7595
sp|Q9KTE0|MIAB_VIBCH tRNA-2-methylthio-N(6)-dimethylallyladenosine synthase OS=Vibrio cholerae serotype O1 (strain ATCC 39315 / El Tor Inaba N16961) OX=243277 GN=miaB PE=3 SV=1: 0.7595
sp|P85808|TRP7_RHOPR Tachykinin-related peptide 7 OS=Rhodnius prolixus OX=13249 PE=1 SV=1: 0.7575
sp|Q80UY1|CARME_MOUSE Carnosine N-methyltransferase OS=Mus musculus OX=10090 GN=Carnmt1 PE=1 SV=1: 0.7549
sp|P42997|OXYF_SCYCA Phasvatocin OS=Scyliorhinus canicula OX=7830 PE=1 SV=1: 0.7454
sp|P42996|OXYA_SCYCA Asvatocin OS=Scyliorhinus canicula OX=7830 PE=1 SV=1: 0.7390
sp|P69128|OXYT_CYPCA Vasotocin OS=Cyprinus carpio OX=7962 PE=1 SV=1: 0.7388
sp|P69129|OXYT_PETMA Vasotocin OS=Petromyzon marinus OX=7757 PE=1 SV=1: 0.7388
sp|O67019|PDXA_AQUAE 4-hydroxythreonine-4-phosphate dehydrogenase OS=Aquifex aeolicus (strain VF5) OX=224324 GN=pdxA PE=3 SV=1: 0.7383
sp|Q5V5M5|TMCA_HALMA tRNA(Met) cytidine acetyltransferase TmcA OS=Haloarcula marismortui (strain ATCC 43049 / DSM 3752 / JCM 8966 / VKM B-1809) OX=272569 GN=tmcA PE=3 SV=1: 0.7354
sp|Q9VTU1|APG12_DROME Autophagy protein 12-like OS=Drosophila melanogaster OX=7227 GN=Atg12 PE=3 SV=3: 0.7342
sp|P43000|OXYV_SQUAC Valitocin OS=Squalus acanthias OX=7797 PE=1 SV=1: 0.7338
sp|E7FGY2|GID8B_DANRE Glucose-induced degradation protein 8-B homolog OS=Danio rerio OX=7955 GN=gid8b PE=2 SV=2: 0.7322
sp|P19118|TPIS_NICPL Triosephosphate isomerase, cytosolic (Fragment) OS=Nicotiana plumbaginifolia OX=4092 PE=1 SV=1: 0.7308
sp|B1I8B3|CH10_STRPI Co-chaperonin GroES OS=Streptococcus pneumoniae (strain Hungary19A-6) OX=487214 GN=groES PE=3 SV=1: 0.7269
sp|Q97NV3|CH10_STRPN Co-chaperonin GroES OS=Streptococcus pneumoniae serotype 4 (strain ATCC BAA-334 / TIGR4) OX=170187 GN=groES PE=1 SV=1: 0.7269
sp|C1CML8|CH10_STRZP Co-chaperonin GroES OS=Streptococcus pneumoniae (strain P1031) OX=488223 GN=groES PE=3 SV=1: 0.7269
sp|P42994|OXYT_RAJCL Glumitocin OS=Raja clavata OX=7781 PE=1 SV=1: 0.7245
sp|A5EX69|RPOC_DICNV DNA-directed RNA polymerase subunit beta' OS=Dichelobacter nodosus (strain VCS1703A) OX=246195 GN=rpoC PE=3 SV=1: 0.7244
sp|P68125|CCKN_DASVI Cholecystokinin OS=Dasyurus viverrinus OX=9279 GN=CCK PE=1 SV=1: 0.7220
sp|P0C229|CCKN_MACGI Cholecystokinin OS=Macropus giganteus OX=9317 GN=CCK PE=1 SV=1: 0.7220
sp|P68126|CCKN_NOTEU Cholecystokinin OS=Notamacropus eugenii OX=9315 GN=CCK PE=1 SV=1: 0.7220
sp|Q7VKI8|MSRP_HAEDU Protein-methionine-sulfoxide reductase catalytic subunit MsrP OS=Haemophilus ducreyi (strain 35000HP / ATCC 700724) OX=233412 GN=msrP PE=3 SV=1: 0.7214
sp|P81814|ALL11_CARMA Carcinustatin-11 OS=Carcinus maenas OX=6759 PE=1 SV=1: 0.7203
sp|O67373|MURC_AQUAE UDP-N-acetylmuramate--L-alanine ligase OS=Aquifex aeolicus (strain VF5) OX=224324 GN=murC PE=3 SV=1: 0.7191
sp|Q862Z2|ASB5_RABIT Ankyrin repeat and SOCS box protein 5 OS=Oryctolagus cuniculus OX=9986 GN=ASB5 PE=2 SV=1: 0.7187
sp|B5FKD6|MALT_SALDC HTH-type transcriptional regulator MalT OS=Salmonella dublin (strain CT_02021853) OX=439851 GN=malT PE=3 SV=1: 0.7186
sp|P38642|UF04_MOUSE Unknown protein from 2D-PAGE of fibroblasts (Fragment) OS=Mus musculus OX=10090 PE=1 SV=1: 0.7176
sp|A0QCX7|ATPG_MYCA1 ATP synthase gamma chain OS=Mycobacterium avium (strain 104) OX=243243 GN=atpG PE=3 SV=1: 0.7167
sp|Q73X58|ATPG_MYCPA ATP synthase gamma chain OS=Mycolicibacterium paratuberculosis (strain ATCC BAA-968 / K-10) OX=262316 GN=atpG PE=3 SV=1: 0.7167
sp|P41860|FAR5_CALVO CalliFMRFamide-5 OS=Calliphora vomitoria OX=27454 PE=1 SV=1: 0.7167
sp|B3EWK2|FAR6_DELRA FMRFamide-like neuropeptide APGQDFMRF-amide OS=Delia radicum OX=30064 PE=1 SV=1: 0.7167
sp|Q7ZWD4|UAP1L_DANRE UDP-N-acetylhexosamine pyrophosphorylase-like protein 1 OS=Danio rerio OX=7955 GN=uap1l1 PE=2 SV=1: 0.7166
sp|P86980|ALDOA_GADMO Fructose-bisphosphate aldolase A (Fragment) OS=Gadus morhua OX=8049 GN=ALDOA PE=1 SV=1: 0.7105
sp|Q9LZB2|GDL74_ARATH GDSL esterase/lipase At5g03980 OS=Arabidopsis thaliana OX=3702 GN=At5g03980 PE=2 SV=1: 0.7067
sp|P84355|PVK2_MUSDO Periviscerokinin-2 OS=Musca domestica OX=7370 PE=1 SV=1: 0.7053
sp|P42999|OXYA_SQUAC Aspartocin OS=Squalus acanthias OX=7797 PE=1 SV=1: 0.7049
sp|Q04SA1|UPPP_LEPBJ Undecaprenyl-diphosphatase OS=Leptospira borgpetersenii serovar Hardjo-bovis (strain JB197) OX=355277 GN=uppP PE=3 SV=1: 0.7045
sp|Q050C7|UPPP_LEPBL Undecaprenyl-diphosphatase OS=Leptospira borgpetersenii serovar Hardjo-bovis (strain L550) OX=355276 GN=uppP PE=3 SV=1: 0.7045
sp|P54713|NDUAA_CANLF NADH dehydrogenase [ubiquinone] 1 alpha subcomplex subunit 10 (Fragment) OS=Canis lupus familiaris OX=9615 GN=NDUFA10 PE=1 SV=1: 0.7043
sp|B3EWP5|BCNGR_LACPA Bacteriocin BacGR (Fragment) OS=Lacticaseibacillus paracasei OX=1597 PE=1 SV=1: 0.7039
sp|Q92EG7|AROE_LISIN Shikimate dehydrogenase (NADP(+)) OS=Listeria innocua serovar 6a (strain ATCC BAA-680 / CLIP 11262) OX=272626 GN=aroE PE=3 SV=1: 0.7035
sp|Q8KRR5|4OT_PSEFL 2-hydroxymuconate tautomerase OS=Pseudomonas fluorescens OX=294 GN=nahJ PE=3 SV=3: 0.7009
sp|P52923|AIF1_YEAST Apoptosis-inducing factor 1 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=AIF1 PE=1 SV=1: 0.6998
sp|P69056|NEU1_BALPH Oxytocin OS=Balaenoptera physalus OX=9770 GN=OXT PE=1 SV=1: 0.6991
sp|P69044|NEU1_HIPAM Oxytocin OS=Hippopotamus amphibius OX=9833 GN=OXT PE=1 SV=1: 0.6991
sp|P69043|NEU1_RABIT Oxytocin OS=Oryctolagus cuniculus OX=9986 GN=OXT PE=1 SV=1: 0.6991
sp|P69057|NEU1_TACAC Oxytocin OS=Tachyglossus aculeatus aculeatus OX=49271 GN=OXT PE=1 SV=1: 0.6991
sp|P69058|OXYT_HYDCO Oxytocin OS=Hydrolagus colliei OX=7873 PE=1 SV=1: 0.6991
sp|Q54CH6|IF4AX_DICDI Putative eukaryotic initiation factor 4A-like protein OS=Dictyostelium discoideum OX=44689 GN=DDB_G0292914 PE=5 SV=1: 0.6949
sp|Q7VGU6|MINE_HELHP Cell division topological specificity factor OS=Helicobacter hepaticus (strain ATCC 51449 / 3B1) OX=235279 GN=minE PE=3 SV=1: 0.6939
sp|B9KZY3|RS19_THERP Small ribosomal subunit protein uS19 OS=Thermomicrobium roseum (strain ATCC 27502 / DSM 5159 / P-2) OX=309801 GN=rpsS PE=3 SV=1: 0.6934
sp|Q81J33|RS17_BACCR Small ribosomal subunit protein uS17 OS=Bacillus cereus (strain ATCC 14579 / DSM 31 / CCUG 7414 / JCM 2152 / NBRC 15305 / NCIMB 9373 / NCTC 2599 / NRRL B-3711) OX=226900 GN=rpsQ PE=3 SV=1: 0.6914
sp|Q97FS5|PYRF_CLOAB Orotidine 5'-phosphate decarboxylase OS=Clostridium acetobutylicum (strain ATCC 824 / DSM 792 / JCM 1419 / IAM 19013 / LMG 5710 / NBRC 13948 / NRRL B-527 / VKM B-1787 / 2291 / W) OX=272562 GN=pyrF PE=3 SV=1: 0.6914
sp|Q64749|PKG1_ADEG1 Packaging protein 1 OS=Fowl adenovirus A serotype 1 (strain CELO / Phelps) OX=10553 GN=IVa2 PE=3 SV=1: 0.6896
sp|B3A0E2|FAR4_AUSGA Extended FMRFamide-4 OS=Austrophasma gansbaaiense OX=253136 PE=1 SV=1: 0.6895
sp|B3A0A3|FAR4_AUSRA Extended FMRFamide-4 OS=Austrophasma rawsonvillense OX=253137 PE=1 SV=1: 0.6895
sp|B3A0C3|FAR4_HEMMO Extended FMRFamide-4 OS=Hemilobophasma montaguense OX=253130 PE=1 SV=1: 0.6895
sp|B3A065|FAR4_KARBI Extended FMRFamide-4 OS=Karoophasma biedouwense OX=253133 PE=1 SV=1: 0.6895
sp|B0M8U3|FAR4_KARBO Extended FMRFamide-4 OS=Karoophasma botterkloofense OX=253132 PE=1 SV=1: 0.6895
sp|B3A084|FAR4_LOBRE Extended FMRFamide-4 OS=Lobatophasma redelinghuysense OX=253128 PE=1 SV=1: 0.6895
sp|B0M3C9|FAR4_MANKU Extended FMRFamide-4 OS=Mantophasma kudubergense OX=1037657 PE=1 SV=1: 0.6895
sp|B0M2T4|FAR4_NAMOO Extended FMRFamide-4 OS=Namaquaphasma ookiepense OX=409167 PE=1 SV=1: 0.6895
sp|B3A0J9|FAR4_PACBA Extended FMRFamide-4 OS=Pachyphasma brandbergense OX=1041430 PE=1 SV=1: 0.6895
sp|B3A0G1|FAR4_PRAMA Extended FMRFamide-4 OS=Praedatophasma maraisi OX=409170 PE=1 SV=1: 0.6895
sp|B0M3A7|FAR4_STRNA Extended FMRFamide-4 OS=Striatophasma naukluftense OX=1041429 PE=1 SV=1: 0.6895
sp|O04130|SERA2_ARATH D-3-phosphoglycerate dehydrogenase 2, chloroplastic OS=Arabidopsis thaliana OX=3702 GN=PGDH2 PE=1 SV=2: 0.6881
sp|P85458|FAR14_LUCCU FMRFamide-14 OS=Lucilia cuprina OX=7375 PE=1 SV=1: 0.6876
sp|A2RUV9|AEBP1_RAT Adipocyte enhancer-binding protein 1 OS=Rattus norvegicus OX=10116 GN=Aebp1 PE=2 SV=1: 0.6873
sp|P81740|TRP8_RHYMA Tachykinin-related peptide 8 OS=Rhyparobia maderae OX=36963 PE=1 SV=1: 0.6872
sp|P77746|YBDO_ECOLI Uncharacterized HTH-type transcriptional regulator YbdO OS=Escherichia coli (strain K12) OX=83333 GN=ybdO PE=2 SV=1: 0.6867
sp|P41861|FAR6_CALVO CalliFMRFamide-6 OS=Calliphora vomitoria OX=27454 PE=1 SV=1: 0.6865
sp|P0DJK2|BPPAD_BOTCO Bradykinin-potentiating peptide 10d OS=Bothrops cotiara OX=8727 PE=1 SV=1: 0.6860
sp|P0DJK5|BPPAD_BOTFO Bradykinin-potentiating peptide 10d OS=Bothrops fonsecai OX=157549 PE=1 SV=1: 0.6860
sp|A9VJW2|TRPB_BACMK Tryptophan synthase beta chain OS=Bacillus mycoides (strain KBAB4) OX=315730 GN=trpB PE=3 SV=1: 0.6858
sp|L0E2U6|PHQG_PENFE NmrA-like family domain-containing oxidoreductase phqG OS=Penicillium fellutanum OX=70095 GN=phqG PE=3 SV=1: 0.6840
sp|Q3BMT4|HEM6_XANE5 Oxygen-dependent coproporphyrinogen-III oxidase OS=Xanthomonas euvesicatoria pv. vesicatoria (strain 85-10) OX=316273 GN=hemF PE=3 SV=1: 0.6833
sp|P0DL04|BPPBI_BOTJA Bradykinin-potentiating peptide 11i OS=Bothrops jararaca OX=8724 PE=1 SV=1: 0.6833
sp|C4R8S9|NOP9_KOMPG Nucleolar protein 9 OS=Komagataella phaffii (strain GS115 / ATCC 20864) OX=644223 GN=NOP9 PE=3 SV=1: 0.6798
sp|P85067|OJP1_SEPOF Ovarian jelly-peptide 1 OS=Sepia officinalis OX=6610 PE=1 SV=1: 0.6788
sp|P84897|BRK6_PITHY [Val1,Hyp2,Thr6]-bradykinyl-Gln OS=Pithecopus hypochondrialis OX=317381 PE=1 SV=1: 0.6784
sp|A1SQS4|RS6_NOCSJ Small ribosomal subunit protein bS6 OS=Nocardioides sp. (strain ATCC BAA-499 / JS614) OX=196162 GN=rpsF PE=3 SV=1: 0.6783
sp|P85466|FAR9_SARBU FMRFamide-9 OS=Sarcophaga bullata OX=7385 PE=1 SV=1: 0.6776
sp|F4IKM1|TRNHB_ARATH Tropinone reductase homolog At2g29340 OS=Arabidopsis thaliana OX=3702 GN=At2g29340 PE=2 SV=1: 0.6772
sp|Q28TS9|RUVB_JANSC Holliday junction branch migration complex subunit RuvB OS=Jannaschia sp. (strain CCS1) OX=290400 GN=ruvB PE=3 SV=1: 0.6769
```

## Как запустить?

- Установить зависимости из requirenments.txt
- Преобразовать файл с помощью read.py (сейчас работает только с архивом uniprot_sprot.fasta.gz)
- Запустить интерпретатор с помощью python3 main.py (ускорено только с numpy)
- Запустить интерпретатор с помощью python3 main_numba.py (ускорено только с numpy, numba etc., работает в несколько десятков раз быстрее)

Пример запуска из командной строки
```
(.venv) PS C:\programming\bioinformatics> python main_numba.py "MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLDAKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHLEKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWKFTPL"
```

## P.S.

В ветке [ark](https://github.com/ark2016/bioinformatics/tree/ark) лежит один из подходов, онсованный на B-spline, но он оказался намного дольше по времени, хотя теоретически может быть лучше для сравнения двух строк, за счёт своей природы.



