scenes = [#"scenes/01_simple_scenes/bunny_cornell.xml",
          "scenes/01_simple_scenes/instanced_objects.xml", #(!)
          "scenes/03_dragon/dragon.xml",                   #(!)
          "scenes/04_bonsai/bonsai.xml",                   #(!)
          "scenes/05_sponza/sponza_c.xml",                 #(!)
          "scenes/06_cry_sponza/cry_sponza_c.xml",         #(!)
          "scenes/06_cry_sponza2/cry_sponza.xml",          #(!)
          #"scenes/07_bistro/bistro_c.xml",
          #"scenes/08_bath/bathroom_c.xml",
          #"scenes/09_hair_balls/hairballs.xml",
          "scenes/10_audi_a8/audi_a8.xml",
          "scenes/11_audi_a8_opt/audi_a8_opt.xml",         #(!)
          "scenes/12_human_sophi/human_sophi.xml",         #(!)
          "scenes/13_hero/hero.xml",
          "scenes/14_village/village.xml",
          ]



options = [#"-accel Embree -build embree",
           #"-accel BVH4Common -build cbvh_embree2  -layout DepthFirst -gpu 0",
           #"-accel BVH4Half   -build cbvh_embree2  -layout DepthFirst -gpu 0",
           #"-accel BVH2Common -build cbvh_embree2  -layout DepthFirst -gpu 0",
           #"-accel BVH2_STAS  -build cbvh_embree2  -layout DepthFirst -gpu 0",
           #"-accel BVH2_LOFT  -build cbvh_embree2  -layout DepthFirst -gpu 0",
           #"-accel BVH2_LRFT  -build cbvh_embree2  -layout DepthFirst -gpu 0",
           #"-accel BVH2Fat32  -build cbvh_embree2  -layout DepthFirst -gpu 0",
           ##"-accel BVH2Fat32  -build cbvh_embree2  -layout Clusterized31 -gpu 0",
           ##"-accel BVH2Fat32  -build cbvh_embree2  -layout SuperTreelet4 -gpu 0",
           ##"-accel BVH2Fat32  -build cbvh_embree2  -layout SuperTreeletAlignedMerged4 -gpu 0",
           #"-accel BVH2Fat16      -build cbvh_embree2  -layout SuperTreelet4 -gpu 0",
           #"-accel BVH2FatCompact -build cbvh_embree2  -layout SuperTreeletAlignedMerged4 -gpu 0",
           ##"-accel ShtRT          -build cbvh_embree2  -layout DepthFirst -gpu 0",
           #"-accel BVH2_SHRT64       -build cbvh_embree2  -layout DepthFirst -gpu 0",
           #"-accel BVH4_SHRT         -build cbvh_embree2  -layout DepthFirst -gpu 0",
           #"-accel BVH4_SHRT_HALF    -build cbvh_embree2  -layout DepthFirst -gpu 0",
           #"-accel BVH2StacklessIfIf -build cbvh_embree2 -layout DepthFirst -gpu 0",
           "-accel RTX            -build RTX -gpu 1",       
           "-accel BVH4_32        -build cbvh_embree2 -layout DepthFirst -gpu 1",
           "-accel BVH4_16        -build cbvh_embree2 -layout DepthFirst -gpu 1",
           "-accel BVH4_16_2D     -build cbvh_embree2 -layout DepthFirst -gpu 1",
           "-accel BVH4_16_2DS    -build cbvh_embree2 -layout DepthFirst -gpu 1",
           #"-accel BVH2Common        -build cbvh_embree2 -layout DepthFirst -gpu 1",
           #"-accel BVH2_STAS         -build cbvh_embree2 -layout DepthFirst -gpu 1",
           #"-accel BVH2_LOFT         -build cbvh_embree2 -layout DepthFirst -gpu 1",
           #"-accel BVH2_LRFT         -build cbvh_embree2 -layout DepthFirst -gpu 1",
           #"-accel BVH2_STACKLESS_II -build cbvh_embree2 -layout DepthFirst -gpu 1",
           #"-accel BVH2_STACKLESS_WW -build cbvh_embree2 -layout DepthFirst -gpu 1", 
           "-accel BVH2Fat32         -build cbvh_embree2 -layout DepthFirst -gpu 1",
           "-accel BVH2Fat16         -build cbvh_embree2 -layout DepthFirst -gpu 1",
           "-accel BVH2Fat16TB       -build cbvh_embree2 -layout DepthFirst -gpu 1",
           ##"-accel BVH2Fat16         -build cbvh_embree2 -layout TreeletBased2 -gpu 1",
           ##"-accel BVH2Fat16         -build cbvh_embree2 -layout TreeletBased4 -gpu 1",
           ##"-accel BVH2Fat16         -build cbvh_embree2 -layout Clusterized31 -gpu 1",
           ##"-accel BVH2Fat06 -build cbvh_embree2 -layout SuperTreeletAlignedMerged8 -gpu 1",
           #"-accel BVH2Fat06 -build cbvh_embree2 -layout SuperTreeletAlignedMerged4 -gpu 1",
           ##"-accel BVH2Fat06 -build cbvh_embree2 -layout SuperTreeletAlignedMerged2 -gpu 1",
           #"-accel BVH2_SHRT32 -build cbvh_embree2  -layout DepthFirst -gpu 1",
           #"-accel BVH2_SHRT64 -build cbvh_embree2  -layout DepthFirst -gpu 1",
           #"-accel BVH4_SHRT   -build cbvh_embree2  -layout DepthFirst -gpu 1",
           #"-accel BVH4_SHRT_HALF -build cbvh_embree2 -layout DepthFirst -gpu 1",
           #"-accel RTX            -build RTX -gpu 1",
           ]

#options = ["-accel RTX       -build RTX          -layout DepthFirst -gpu 1",
#           "-accel BVH2Fat16 -build nanort2      -layout DepthFirst -gpu 1",
#           "-accel BVH2Fat16 -build nanort4      -layout DepthFirst -gpu 1",
#           "-accel BVH2Fat16 -build cbvh_lbvh    -layout DepthFirst -gpu 1",
#           "-accel BVH2Fat16 -build cbvh_hq4     -layout DepthFirst -gpu 1",
#           "-accel BVH2Fat16 -build cbvh_med4    -layout DepthFirst -gpu 1",
#           "-accel BVH2Fat16 -build cbvh_embree2 -layout DepthFirst -gpu 1",
#           "-accel BVH2Fat16 -build cbvh_embree4 -layout DepthFirst -gpu 1",
#           "-accel BVH2Fat06 -build cbvh_hq4     -layout SuperTreeletAlignedMerged4 -gpu 1",
#           "-accel BVH2Fat06 -build cbvh_med4    -layout SuperTreeletAlignedMerged4 -gpu 1",
#           "-accel BVH2Fat06 -build cbvh_embree2 -layout SuperTreeletAlignedMerged4 -gpu 1",
#           "-accel BVH2Fat06 -build cbvh_embree4 -layout SuperTreeletAlignedMerged4 -gpu 1",
#           "-accel BVH4_16   -build cbvh_hq4     -layout DepthFirst -gpu 1",
#           "-accel BVH4_16   -build cbvh_med4    -layout DepthFirst -gpu 1",
#           "-accel BVH4_16   -build cbvh_embree2 -layout DepthFirst -gpu 1",
#           "-accel BVH4_16   -build cbvh_embree4 -layout DepthFirst -gpu 1",
#           ]

#options = ["-render RT -accel RTX            -build RTX -gpu 1",
#           "-render RT -accel BVH4Common     -build cbvh_embree2 -layout DepthFirst -gpu 1", 
#           "-render RT -accel BVH4Half       -build cbvh_embree2 -layout DepthFirst -gpu 1",
#           "-render RT -accel BVH4_SHRT_HALF -build cbvh_embree2 -layout DepthFirst -gpu 1",
#           "-render RT -accel BVH2_SHRT64    -build cbvh_embree2  -layout DepthFirst -gpu 1",
#           "-render RT -accel BVH2_LOFT      -build cbvh_embree2 -layout DepthFirst -gpu 1",
#           "-render RT -accel BVH2Fat32      -build cbvh_embree2 -layout DepthFirst -gpu 1",
#           "-render RT -accel BVH2Fat16      -build cbvh_embree2 -layout DepthFirst -gpu 1",
#           "-render RT -accel BVH2Fat06      -build cbvh_embree2 -layout SuperTreeletAlignedMerged4 -gpu 1", 
#           ]   


options = ["-render AO -accel RTX            -build RTX -gpu 1",
           "-render AO -accel BVH4_32        -build cbvh_embree2 -layout DepthFirst -gpu 1", 
           "-render AO -accel BVH4_16        -build cbvh_embree2 -layout DepthFirst -gpu 1",
           "-render AO -accel BVH4_SHRT_HALF -build cbvh_embree2 -layout DepthFirst -gpu 1",
           "-render AO -accel BVH2_SHRT64    -build cbvh_embree2 -layout DepthFirst -gpu 1",
           "-render AO -accel BVH2_STACKLESS -build cbvh_embree2 -layout DepthFirst -gpu 1", # enable for AO?
           "-render AO -accel BVH2_STAS      -build cbvh_embree2 -layout DepthFirst -gpu 1", # enable for AO?
           "-render AO -accel BVH2_LOFT      -build cbvh_embree2 -layout DepthFirst -gpu 1",
           "-render AO -accel BVH2Fat32      -build cbvh_embree2 -layout DepthFirst -gpu 1",
           "-render AO -accel BVH2Fat16      -build cbvh_embree2 -layout DepthFirst -gpu 1",
           "-render AO -accel BVH2Fat06      -build cbvh_embree2 -layout SuperTreeletAlignedMerged4 -gpu 1",
           ]  

#options = ["-render RT -accel RTX -build RTX -gpu 1"]                 

#options = ["-render EyeRay -accel RTX       -build RTX -gpu 1",
#           "-render EyeRay -accel BVH2_LRFT -build cbvh_lbvh -layout DepthFirst -gpu 1",
#           "-render EyeRay -accel BVH2_LOFT -build cbvh_lbvh -layout DepthFirst -gpu 1"
#           ] 

#options = ["-render EyeRay -accel Embree    -build embree    -layout DepthFirs  -gpu 0",
#           "-render EyeRay -accel BVH2_LRFT -build cbvh_lbvh -layout DepthFirst -gpu 0",
#           "-render EyeRay -accel BVH2_LOFT -build cbvh_lbvh -layout DepthFirst -gpu 0"
#           ] 

#scenes = ["scenes/03_dragon/dragon2.xml",                
#          "scenes/06_cry_sponza/cry_sponza_c2.xml",              
#          "scenes/16_conference_collapsed/conference_c2.xml",
#          "scenes/16_conference_collapsed/conference_c3.xml",
#          ]

#options = ["-accel BVH2Fat32 -build cbvh_embree2 -layout DepthFirst -gpu 0"]


#options = ["-render EyeRay -accel BVH2_LOFT -build cbvh_embree4 -gpu 0",
#           "-render EyeRay -accel BVH2_LOFT -build cbvh_hq4     -gpu 0",
#           "-render EyeRay -accel BVH2_LOFT -build cbvh_lbvh    -gpu 0"] 

#options = ["-render EyeRay -accel BVH2_LOFT -build cbvh_embree4 -gpu 0",
#           "-render EyeRay -accel BVH2_LRFT -build cbvh_lbvh    -gpu 1",
#           "-render EyeRay -accel BVH2_LOFT -build cbvh_lbvh    -gpu 1",
#           ] 

options = ["-render EyeRay -accel BVH2Fat32 -build nanort2            -gpu 1",
           "-render EyeRay -accel BVH2Fat32 -build cbvh_embree2       -gpu 1",
           "-render EyeRay -accel BVH2Fat32 -build cbvh_embree_fast2  -gpu 1",
           "-render EyeRay -accel BVH2Fat32 -build cbvh_hq2           -gpu 1",
           "-render EyeRay -accel BVH2Fat32 -build cbvh_med4          -gpu 1",
           "-render EyeRay -accel BVH2Fat32 -build cbvh_lbvh          -gpu 1",
           "-render EyeRay -accel BVH2Fat32 -build cbvh_lbvh_gpu      -gpu 1",
           "-render EyeRay -accel BVH2_LRFT -build cbvh_lbvh_gpu      -gpu 1"] 

timingsFile  = open("run_tests.sh", "w")
timingsFile2 = open("run_tests_bf.sh", "w")
timingsFile.write("mkdir test_images\n")
for (scene, scnId) in zip(scenes, range(len(scenes))):
  (A,B,C) = scene.split("/")
  #print(A,B,C)
  image = B if B != "01_simple_scenes" else "0" + str(scnId) + "_" + C.split(".")[0]
  for (option,opid) in zip(options, range(len(options))):   
    timingsFile.write("./render_app -scene \"" + scene + "\"" + " -out test_images/x" + image + f'{opid:02d}' + ".bmp " + option + " -test 1 -launch_id " + str(opid) + " -width 2048 -height 2048 \n")
    timingsFile2.write("./render_app -render bfrt -scene \"" + scene + "\"" + " -out test_images/x" + image + f'{opid:02d}' + ".bmp " + option + " -test 1 -launch_id " + str(opid) + " -width 2048 -height 2048 \n")

timingsFile.write ("./render_app -scene \"scenes/01_simple_scenes/instanced_objects.xml\" -out test_images/x99_instanced_objects08.bmp -accel Overhead -build cbvh_embree2 -gpu 1 -test 1 -launch_id 99 -width 2048 -height 2048") 
timingsFile2.write("./render_app -scene \"scenes/01_simple_scenes/instanced_objects.xml\" -out test_images/x99_instanced_objects08.bmp -accel Overhead -build cbvh_embree2 -gpu 1 -test 1 -launch_id 99 -render bfrt -width 2048 -height 2048") 

timingsFile.close()
timingsFile2.close()
