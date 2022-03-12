
# import some common libraries
import numpy as np
import cv2
import matplotlib.pyplot as plt
from detectron2 import model_zoo
from detectron2.engine import DefaultPredictor
from detectron2.config import get_cfg
from detectron2.utils.visualizer import Visualizer
from detectron2.data import MetadataCatalog
from detectron2.utils.visualizer import ColorMode



if __name__ == "__main__":
    cfg = get_cfg()
    cfg.merge_from_file("./configs/COCO-InstanceSegmentation/mask_rcnn_R_50_FPN_3x.yaml")
    cfg.MODEL.ROI_HEADS.SCORE_THRESH_TEST = 0.5  #模型阈值
    cfg.MODEL.WEIGHTS = "model_final.pth"
    cfg.MODEL.ROI_HEADS.NUM_CLASSES = 1

    predictor = DefaultPredictor(cfg)
    im = cv2.imread('1.png')
    outputs = predictor(im)

    # 实例分割结果显示
    v = Visualizer(im[:, :, ::-1], MetadataCatalog.get(cfg.DATASETS.TRAIN[0]), scale=1.2)
    v = v.draw_instance_predictions(outputs["instances"].to("cpu"))
    plt.figure()
    plt.imshow(v.get_image()[:, :, ::-1])


    # 二值化结果显示
    height=im.shape[0]   
    width=im.shape[1] 
    img_seg = np.zeros((height,width))
    img_seg = img_seg.astype(bool)
    if len(outputs["instances"])>0:
        instance = outputs["instances"].to("cpu").pred_masks[0,:,:].detach().numpy()
        for m in range(1,len(outputs["instances"])):
            instance = np.bitwise_or(instance,outputs["instances"].to("cpu").pred_masks[m,:,:].detach().numpy())
        img_seg = np.bitwise_or(instance, img_seg)
    
    plt.figure()
    plt.imshow(img_seg,cmap="gray")
    plt.show()
        

