scenes = [#"scenes/01_simple_scenes/bunny_cornell.xml",
          #"scenes/01_simple_scenes/instanced_objects.xml",
          "scenes/03_dragon/dragon.xml",
          #"scenes/04_bonsai/bonsai.xml",
          "scenes/05_sponza/sponza_c.xml",
          "scenes/06_cry_sponza/cry_sponza_c.xml",
          "scenes/06_cry_sponza2/cry_sponza.xml",
          "scenes/07_bistro/bistro_c.xml",
          #"scenes/08_bath/bathroom_c.xml",
          #"scenes/09_hair_balls/hairballs.xml",

          #"scenes/10_audi_a8/audi_a8.xml",
          #"scenes/11_audi_a8_opt/audi_a8_opt.xml",
          #"scenes/12_human_sophi/human_sophi.xml",
          #"scenes/13_hero/hero.xml",
          #"scenes/14_village/village.xml",
          ]


images = [#"01_bunny",
          #"02_inst_obj",
          "03_dragon",
          #"04_bonsai",
          "05_sponza",
          "06_cry_sponza",
          "06_cry_sponza2",
          "07_bistro",
          #"08_bath",
          #"09_hairballs",
          
          #"10_audi_a8",
          #"11_audi_a8_opt",
          #"12_human_sophi",
          #"13_human_hero",
          #"14_village"
          ]


#scenes = scenes[3:4]
#images = images[3:4]

optionsTime = [#"-accel Embree     -build embree",
               #"-accel BVH4Common -build cbvh_embree2  -layout DepthFirst -gpu 0",
               #"-accel BVH4Half   -build cbvh_embree2  -layout DepthFirst -gpu 0",
               #"-accel BVH2Common -build cbvh_embree2  -layout DepthFirst -gpu 0",
               #"-accel BVH2_LOFT  -build cbvh_embree2  -layout DepthFirst -gpu 0",
               #"-accel BVH2_LRFT  -build cbvh_embree2  -layout DepthFirst -gpu 0",
               #"-accel BVH2_STAS  -build cbvh_embree2  -layout DepthFirst -gpu 0",
               #"-accel BVH2Fat    -build cbvh_embree2  -layout DepthFirst -gpu 0",
               #"-accel BVH2Fat16  -build cbvh_embree2  -layout SuperTreelet4 -gpu 0",
               #"-accel BVH2FatCompact -build cbvh_embree2  -layout SuperTreeletAlignedMerged4 -gpu 0",
               #"-accel ShtRT          -build cbvh_embree2  -layout DepthFirst",
               #"-accel BVH2StacklessIfIf -build cbvh_embree2 -layout DepthFirst -gpu 0",
               ]
               
                              
optionsMetr1 = [
               "-accel BVH4Common -build cbvh_embree2 -layout DepthFirst -gpu 0",
               #"-accel BVH4Half   -build cbvh_embree2 -layout DepthFirst -gpu 0",
               #"-accel BVH2Common -build cbvh_embree2 -layout DepthFirst -gpu 0",
               #"-accel BVH2_LOFT  -build cbvh_embree2 -layout DepthFirst -gpu 0",
               #"-accel BVH2_STAS  -build cbvh_embree2 -layout DepthFirst -gpu 0",
               "-accel BVH2Fat    -build cbvh_embree2 -layout DepthFirst -gpu 0",
               #"-accel BVH2Fat16  -build cbvh_embree2 -layout SuperTreelet4 -gpu 0",
               #"-accel BVH2FatCompact -build cbvh_embree2 -layout SuperTreeletAlignedMerged4 -gpu 0",
               #"-accel ShtRT          -build cbvh_embree2  -layout DepthFirst",
               #"-accel BVH2StacklessIfIf -build cbvh_embree2 -layout DepthFirst -gpu 0",
               ]               

optionsMetr2 = [
               #"-accel BVH2Common -build cbvh_embree2",
               #"-accel BVH2StacklessNaive -build cbvh_embree2",
               #"-accel ShtRT -build cbvh_embree2", 
               #"-accel ShtRT     -build cbvh_embree2 -layout DepthFirst",
               # "-accel NanoRTExt  -build NanoRTExt -layout DepthFirst",
               # "-accel BVH4Common -build cbvh_embree -layout DepthFirst",
               # "-accel BVH4Common -build cbvh_med4 -layout DepthFirst",
               # "-accel BVH2Fat -build cbvh_embree2 -layout BreadthFirst",
               #"-accel BVH2Fat     -build nanort2      -layout DepthFirst",
               #"-accel BVH2Fat     -build cbvh_embree2 -layout DepthFirst",
               #"-accel BVH2Fat -build cbvh_embree2 -layout OrderedDepthFirst",
               #"-accel BVH2Fat     -build cbvh_embree2 -layout Clusterized31",
               #"-accel BVH2Fat     -build cbvh_embree2 -layout SuperTreelet4",
               #"-accel BVH2FatHalf -build nanort2      -layout DepthFirst",
               #"-accel BVH2FatHalf -build cbvh_embree2 -layout Clusterized31",
               #"-accel BVH2Fat -build cbvh_embree2 -layout TreeletBased8",
               #"-accel BVH2Fat -build cbvh_embree2 -layout TreeletBased4", 
               #"-accel BVH2Fat -build cbvh_embree2 -layout SuperTreelet8",
               #"-accel BVH2Fat -build cbvh_embree2 -layout SuperTreeletAlignedMerged4",
               #"-accel BVH2Fat -build cbvh_embree2 -layout SuperTreeletAlignedMerged8",
               #"-accel BVH2FatHalf    -build cbvh_embree2 -layout SuperTreelet4",
               #"-accel BVH2FatCompact -build cbvh_embree2 -layout SuperTreelet4",
               #"-accel BVH2FatCompact -build cbvh_embree2 -layout SuperTreeletAlignedMerged4",
               ]  

timingsFile = open("run_timings.sh", "w")
for (scene,image) in zip(scenes,images):
  for (option,opid) in zip(optionsTime, range(len(optionsTime))):
    timingsFile.write("./render_app -scene \"" + scene + "\"" + " -out test_images2/x" + image + str(opid) + ".bmp " + option + "\n")
timingsFile.close()

test_modes = [ ("run_metrics_eye", "test_images2/x", "eye","1"   , "w_stats_eye", "-width 2048 -height 2048"), 
               ("run_metrics_dif", "test_images2/y", "ao" ,"0.75", "w_stats_dif", "-width 512 -height 512"), 
               ("run_metrics_aos", "test_images2/z", "ao" ,"0.05", "w_stats_aos", "-width 512 -height 512") ]

opt_modes = [optionsMetr1,optionsMetr2]
optNames  = ["b","l"]
#opt_modes = [optionsMetr2]
#optNames  = ["l"]

for (optionsMetr, opName) in zip(opt_modes,optNames):
  for (shName,imgPrefix,renderType,rayLength,statsFile, whOpt) in test_modes:
    metricsFile = open(shName + "_" + opName + ".sh", "w")
    for (scene,image) in zip(scenes,images):
      for (option,opid) in zip(optionsMetr, range(len(optionsMetr))):
          metricsFile.write("./render_app -scene \""  + scene + "\"" + \
                            " -out "      + imgPrefix  + image + str(opid) + ".bmp" + \
                            " -render "   + renderType + \
                            " -aoraylen " + rayLength  + \
                            " -stats "    + statsFile  + "_" + opName + ".csv" + \
                            " " + whOpt + " "  + option + "-width 2048 -height 2048 \n")
  metricsFile.close()


metricsFile = open("run_ao_ref.sh", "w")
for (scene,image) in zip(scenes,images):
   metricsFile.write("./render_app -scene \""  + scene + "\"" + " -out test_images2/y" + image + ".bmp -render ao -accel Embree -aoraylen 0.1 -aoraynum 256 -width 1920 -height 1080 \n")
