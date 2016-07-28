etc_json = gen.download_etc_data(experiments[0].id)
print("Time Points: "  + str(etc_json["etc"]["timePoints"]))
print("Gene DPU_G0054708:  " + str(etc_json['etc']['genes']["DPU_G0054708"]))
