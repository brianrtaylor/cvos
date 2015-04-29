function box_centers = bbox_to_box_center(bboxes)
box_centers = [cat(1, bboxes.y), cat(1, bboxes.x)];
end