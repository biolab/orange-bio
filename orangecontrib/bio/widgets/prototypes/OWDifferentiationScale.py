"""
<name>Differentiation Scale</name>
<description></description>
<prototype>1</prototype>
"""

from __future__ import absolute_import

import os, sys
import random
from collections import defaultdict
from operator import itemgetter, add

import numpy

import Orange
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from ... import obiDifscale

class OWDifferentiationScale(OWWidget):
    def __init__(self, parent=None, signalManager=None, title="Differentiation Scale"):
        OWWidget.__init__(self, parent, signalManager, title, wantGraph=True)
        
        self.inputs = [("Gene Expression Samples", Orange.data.Table, self.set_data), ("Additional Expression Samples", Orange.data.Table, self.set_additional_data)]
        self.outputs = [("Selected Time Points", Orange.data.Table), ("Additional Selected Time Points", Orange.data.Table)]
        
        self.selected_time_label = 0
        self.auto_commit = 0
        
        self.loadSettings()
        
        self.selection_changed_flag = False
        
        #####
        # GUI
        #####
        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.info_label = OWGUI.widgetLabel(box, "No data on input")
        self.info_label.setWordWrap(True)
        self.info_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        OWGUI.rubber(self.controlArea)
        
        box = OWGUI.widgetBox(self.controlArea, "Selection")
        
        cb = OWGUI.checkBox(box, self, "auto_commit", "Commit on any change",
                            tooltip="Send updated selections automatically",
                            callback=self.commit_if)
        
        b = OWGUI.button(box, self, "Commit",
                         callback=self.commit,
                         tooltip="Send selections on output signals")
        
        OWGUI.setStopper(self, b, cb, "selection_changed_flag",
                         callback=self.commit)
        
        self.connect(self.graphButton, SIGNAL("pressed()"), self.save_graph)
        
        self.scene = QGraphicsScene()
        self.scene_view = DiffScaleView(self.scene, self.mainArea)
        self.scene_view.setRenderHint(QPainter.Antialiasing)
        self.scene_view.setMinimumWidth(300)
        self.mainArea.layout().addWidget(self.scene_view)
        self.connect(self.scene, SIGNAL("selectionChanged()"), self.on_selection_changed)
        self.connect(self.scene_view, SIGNAL("view_resized(QSize)"), lambda size: self.on_view_resized())
        
        self.data = None
        self.additional_data = None
        self.projections1 = []
        self.projections2 = []
        self.labels1 = []
        self.labels2 = []
        
        self.selected_time_samples = [], []
        
        self.controlArea.setMaximumWidth(300)
        self.resize(600, 480)
        
    def clear(self):
        """ Clear the widget state
        """
        self.projections1 = []
        self.projections2 = []
        self.labels1 = []
        self.labels2 = []
        self.clear_selection()
        self.scene.clear()
        
    def clear_selection(self):
        """ Clear the time point selection.
        """
        self.selected_time_samples = [], []
        
    def set_data(self, data = None):
        """ Set the data for the widget.
        """
        self.clear()
        self.data = data
        
    def set_additional_data(self, data=None):
        """ Set an additional data set.
        """
        self.clear()
        self.additional_data = data
        
    def handleNewSignals(self):
        if self.data is not None:
            self.run_projections()
            self.projection_layout()
            self.update_graph()
            
            info_text = """\
Data with {0} genes 
and {1} samples on input.\n""".format(len(self.data),
                 len(self.data.domain.attributes))
            if self.additional_data is not None:
                info_text += """\
Additional data with {0} genes
and  {1} samples on input.""".format(len(self.additional_data),
                                                    len(self.additional_data.domain.attributes))
            self.info_label.setText(info_text)
        else:
            self.send("Selected Time Points", None)
            self.send("Additional Selected Time Points", None)
            self.info_label.setText("No data on input\n")
            
    def run_projections(self):
        """ Run obiDifscale.get_projections with the current inputs.
        """
        self.error()
#        try:
#            attr_set = list(set(a.attributes['time'] for a in data.domain.attributes))
#            self.time_points = obiDifscale.conv(attr_set, ticks=False)
#        except KeyError, ex:
#            self.error("Could not extract time data")
#            self.clear()
#            return
        
        try:
            (self.projections1, self.labels1,
             self.projections2, self.labels2) = \
                obiDifscale.get_projections(self.data, data2=self.additional_data)
        except Exception, ex:
            self.error("Failed to obtain the projections due to: %r" % ex)
            self.clear()
            return
        
    def projection_layout(self):
        """ Compute the layout for the projections.
        """
        if self.projections1: 
            projections = self.projections1 + self.projections2
            projections = numpy.array(projections)
            
            x_min = numpy.min(projections)
            x_max = numpy.max(projections)
            
            # Scale projections 
            projections = (projections - x_min) / ((x_max - x_min) or 1.0)
            projections = list(projections)
            
            labels = self.labels1 + self.labels2
            
            samples = [(attr, self.data) for attr in self.data.domain.attributes] + \
                      ([(attr, self.additional_data) for attr in self.additional_data.domain.attributes] \
                       if self.additional_data is not None else [])
            
            # TODO: handle samples with the same projection
            # the point_layout should return the proj to sample mapping instead
            proj_to_sample = dict([((label, proj), sample) for label, proj, sample \
                                   in zip(labels, projections, samples)])
            self.proj_to_sample = proj_to_sample
            
            time_points = point_layout(labels, projections)
            self.time_points = time_points
            level_height = 20
            all_points = numpy.array(reduce(add, [p for _, p in time_points], []))
            self.all_points = all_points
            
#            all_points[:, 1] *= -level_height
            self.time_samples = [] # samples for time label (same order as in self.time_points)
            
            point_i = 0
            for label, points, in time_points:
                samples = [] 
                for x, y in points:
                    samples.append(proj_to_sample.get((label, x), None))
                self.time_samples.append((label, samples))
            
    def update_graph(self):
        """ Populate the Graphics Scene with the current projections. 
        """
        scene_size_hint = self.scene_view.viewport().size()
        scene_size_hint = QSizeF(max(scene_size_hint.width() - 50, 100),
                                 scene_size_hint.height())
        self.scene.clear()
        
        if self.projections1:
            level_height = 20
            all_points = self.all_points.copy()
            all_points[:, 0] *= scene_size_hint.width()
            all_points[:, 1] *= -level_height
            
            point_i = 0
            centers = []
            z_value = 0
            for label, samples in self.time_samples:
                # Points
                p1 = all_points[point_i]
                points = all_points[point_i: point_i + len(samples), :]
                for (x, y), sample in zip(points, samples):
                    item = GraphicsTimePoint(QRectF(QPointF(x-3, y-3), QSizeF(6, 6)))
                    item.setBrush(QBrush(Qt.black))
                    item.sample = sample
                    item.setToolTip(sample[0].name if sample else "")
                    item.setZValue(z_value)
                    self.scene.addItem(item)
                    point_i += 1
                p2 = all_points[point_i - 1]
                
                # Line over all points
                line = QGraphicsLineItem(QLineF(*(tuple(p1) + tuple(p2))))
                line.setPen(QPen(Qt.black, 2))
                line.setZValue(z_value - 1)
                self.scene.addItem(line)
                
                # Time label on top of the median
                n_points = len(points)
                if n_points % 2:
                    center = points[n_points / 2]
                else:
                    center = (points[n_points / 2] + points[n_points / 2 + 1]) / 2.0
                centers.append(center)
                x, y = center
                text = QGraphicsSimpleTextItem(label)
                w = text.boundingRect().width()
                text.setPos(x - w / 2.0, y - 17.5)
                self.scene.addItem(text)
            
            self.scene.addLine(QLineF(0.0, 0.0, scene_size_hint.width(), 0.0))
            
            polygon = QPolygonF([QPointF(3.0, 0.0),
                                 QPointF(-2.0, -2.0),
                                 QPointF(0.0, 0.0),
                                 QPointF(-2.0, 2.0),
                                 QPointF(3.0, 0.0)])
            
            arrow = QGraphicsPolygonItem(polygon)
            arrow.setBrush(QBrush(Qt.black))
            arrow.setPos(scene_size_hint.width(), 0.0)
            arrow.scale(2, 2)
            self.scene.addItem(arrow)
            
            title = QGraphicsSimpleTextItem("Development (time)")
            font = self.font()
            font.setPointSize(10)
            title.setFont(font)
            w = title.boundingRect().width()
            title.setPos(scene_size_hint.width() - w, -15)
            self.scene.addItem(title)
            
            rects = []
            ticks = []
            axis_label_items = []
            labels = [(center, label) for center, (label, _) in zip(centers, self.time_samples)]
            labels = sorted(labels, key=lambda (c, l): c[0])
            for center, label in labels:
                x, y = center
                item = QGraphicsSimpleTextItem(label)
                w = item.boundingRect().width()
                item.setPos(x - w / 2.0, 4.0)
                rects.append(item.sceneBoundingRect().normalized())
                ticks.append(QPointF(x - w / 2.0, 4.0))
                axis_label_items.append(item)
            
#            rects = SA_axis_label_layout(ticks, rects, max_time=0.5,
#                                         x_factor=scene_size_hint.width() / 50.0,
#                                         y_factor=10,
#                                         random=random.Random(0))

            rects = greedy_scale_label_layout(ticks, rects, spacing=5)
            
            for (tick, label), rect, item in zip(labels, rects, axis_label_items):
                x, y = tick
                self.scene.addLine(x, -2, x, 2)
                if rect.top() - item.pos().y() > 5:
                    self.scene.addLine(x, 2, rect.center().x(), 14.0)
                if rect.top() - item.pos().y() > 15:
                    self.scene.addLine(rect.center().x(), 14.0, rect.center().x(), rect.top())
#                item.setPos(rect.topLeft())
                
#                text = QGraphicsSimpleTextItem(label)
#            for tick, rect, item in zip(ticks, rects, axis_label_items):
                item.setPos(rect.topLeft())
                self.scene.addItem(item)
#                w = text.boundingRect().width()
#                text.setPos(x - w / 2.0, 4)
                # Need to compute axis label layout.
#                self.scene.addItem(text)

            self.scene.setSceneRect(self.scene.itemsBoundingRect().adjusted(-10, -10, 10, 10))

    def on_view_resized(self):
        self.update_graph()
        
    def on_selection_changed(self):
        try:
            selected = self.scene.selectedItems()
        except RuntimeError:
            return
        
        selected_attrs1 = []
        selected_attrs2  =[]
        for point in selected:
            attr, data = point.sample if point.sample else (None, None)
            if data is self.data:
                selected_attrs1.append(attr)
            elif data is self.additional_data:
                selected_attrs2.append(attr)
                
        self.selected_time_samples = selected_attrs1, selected_attrs2
        print self.selected_time_samples
        self.commit_if()
            
    def commit_if(self):
        if self.auto_commit:
            self.commit()
        else:
            self.selection_changed_flag = True
    
    def commit(self):
        if self.data is not None:
            selected1, selected2 = self.selected_time_samples
            attrs1 = [a for a in self.data.domain.attributes \
                      if a in selected1]
            domain = Orange.data.Domain(attrs1, self.data.domain.class_var)
            domain.add_metas(self.data.domain.get_metas())
            data = Orange.data.Table(domain, self.data)
            self.send("Selected Time Points", data)
            
            if self.additional_data is not None:
                attrs2 = [a for a in self.additional_data.domain.attributes \
                          if a in selected2]
                domain = Orange.data.Domain(attrs2, self.additional_data.domain.class_var)
                domain.add_metas(self.additional_data.domain.get_metas())
                data = Orange.data.Table(domain, self.additional_data)
                self.send("Additional Selected Time Points", data)
        else:
            self.send("Selected Time Points", None)
            self.send("Additional Selected Time Points", None)
        self.selection_changed_flag = False
        
    def save_graph(self):
        from Orange.OrangeWidgets.OWDlgs import OWChooseImageSizeDlg
        dlg = OWChooseImageSizeDlg(self.scene, parent=self)
        dlg.exec_()
    
    
class GraphicsTimePoint(QGraphicsEllipseItem):
    def __init__(self, *args):
        QGraphicsEllipseItem.__init__(self, *args)
        self.setFlags(QGraphicsItem.ItemIsSelectable)
        self.setAcceptsHoverEvents(True)
        self._is_hovering = False
        
    def paint(self, painter, option, widget=0):
        if self.isSelected():
            brush = QBrush(Qt.red)
            pen = QPen(Qt.red, 1)
        else:
            brush = QBrush(Qt.darkGray)
            pen = QPen(Qt.black, 1)
        if self._is_hovering:
            brush = QBrush(brush.color().darker(200))
        painter.save()
        painter.setBrush(brush)
        painter.setPen(pen)
        painter.drawEllipse(self.rect())
        painter.restore()
        
    def hoverEnterEvent(self, event):
        self._is_hovering = True
        self.update()
        return QGraphicsEllipseItem.hoverEnterEvent(self, event)
    
    def hoverLeaveEvent(self, event):
        self._is_hovering = False
        self.update()
        return QGraphicsEllipseItem.hoverLeaveEvent(self, event)
        
    
class DiffScaleView(QGraphicsView):
    def resizeEvent(self, event):
        QGraphicsView.resizeEvent(self, event)
        self.emit(SIGNAL("view_resized(QSize)"), event.size())
        

def point_layout(labels, points, label_size_hints=None):
    groups = defaultdict(list)
    for label, point in zip(labels, points):
        groups[label].append(point)
        
    for label, points in list(groups.items()):
        points = sorted(points)
        # TODO: Use label_size_hints for min, max
        groups[label] = (points, (points[0], points[-1]))
    
    sorted_groups = sorted(groups.items(), key=itemgetter(1), reverse=True)
    levels = {}
    curr_level = 1
    label_levels = {}
    while sorted_groups:
        label, (points, (x_min, x_max)) = sorted_groups.pop(-1)
        max_level_pos = levels.get(curr_level, x_min)
        if x_min < max_level_pos:
            curr_level += 1
            sorted_groups.append((label, (points, (x_min, x_max))))
        else:
            label_levels[label] = curr_level
            levels[curr_level] = x_max
            curr_level = 1
            
    for label, (points, _) in list(groups.items()):
        level = float(label_levels[label])
        groups[label] = [(x, level) for x in points]
        
    return list(groups.items())

    
def greedy_scale_label_layout(ticks, rects, spacing=3):
    """ Layout the labels at ticks on a linear scale, by raising the
    overlapping labels.
    
    """
    def adjust_interval(start, end, min_v, max_v):
        """ Adjust (start, end) interval to fit inside the (min_v, max_v).
        """
        if start < min_v:
            return (min_v, min_v + (end - start))
        elif max_v > end:
            return (max_v - (end - start), max_v)
        else:
            return (start, end)
        
    def center_interval(start, end, center):
        """ Center the interval on `center`
        """
        span = end - start
        return centered(center, span)
    
    def centered(center, span):
        """ Return an centered interval with span.
        """
        return (center - span / 2.0, center + span / 2.0)
    
    def contains((start, end), (start1, end1)):
        return start <= start1  and end >= end1
    
    def fit(work, ticks, min_x, max_x):
        """ Fit the work set between min_x and max_x  and centered on the
        ticks, if possible.
        """
        fits = False
        work_set = map(QRectF, work)
        tick_center = sum([r.center().x() for r in work_set]) / len(work_set)
        if len(work_set) == 1:
            if work_set[0].left() >= min_x and work_set[0].right() <= max_x:
                return work_set
            else:
                return []
        
        elif len(work_set) == 2: # TODO: MErge this with the > 2
            w_sum = sum([r.width() for r in work_set]) + spacing
            if w_sum < max_x - min_x:
                r1, r2 = work_set
                interval = centered(tick_center, w_sum)
                
                if not contains((min_x, max_x), interval):
                    interval = adjust_interval(*(interval + (min_x, max_x)))
                    
                if contains((min_x, max_x), interval):
                    r1.moveLeft(interval[0])
                    r2.moveLeft(interval[1] - r2.width())
                    r1.moveTop(r1.top() + 10)
                    r2.moveTop(r2.top() + 10)
                    return work_set
                else:
                    return []
            else:
                return []
        
        elif len(work_set) > 2:
            center = (work_set[0].center().x() + work_set[-1].center().x()) / 2.0
            w_sum = work_set[0].width() / 2.0 + work_set[-1].width() / 2.0 + spacing
            for i, r in enumerate(work_set[1:-1]):
                w_sum += r.width() + spacing
            interval = centered(center, w_sum)
            
            if not contains((min_x, max_x), interval):
                interval = adjust_interval(*(interval + (min_x, max_x)))
                
            if contains((min_x, max_x), interval):
                istart, iend = interval
                rstart, rend = work_set[0], work_set[-1]
                rstart.moveLeft(istart)
                rstart.moveTop(rstart.top() + 10)
                rend.moveLeft(iend - rend.width())
                rend.moveTop(rend.top() + 10)
                istart += rstart.width() / 2.0
                iend -= rend.width() / 2.0
                for r in work_set[1: -1]:
                    r.moveLeft(istart)
                    r.moveTop(r.top() + 20)
                    istart += r.width() + spacing
                return work_set
            else:
                return []
            
    queue = sorted(zip(ticks, rects),
                   key=lambda (t, _): t.x(),
                   reverse=True)
    done = False
    rects = []
    
    min_x = -1e30
    max_x = 1e30
    
    while queue:
        work_set = [queue.pop(-1)]
        set_fits = False
        max_x = queue[-1][1].left() if queue else 1e30
        while not set_fits:
            new_rects = fit(map(itemgetter(1), work_set),
                            map(itemgetter(0), work_set),
                            min_x, max_x)
            if new_rects: # Can the work set be fit.
                set_fits = True
                rects.extend(new_rects)
                min_x = work_set[-1][1].right()
                
            else:
                # Extend the work set with one more label rect
                work_set.append(queue.pop(-1))
                max_x = queue[-1][1].left() if queue else 1e30
    return rects
        
    
if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWDifferentiationScale()
    data = Orange.data.Table(os.path.expanduser("~/Documents/GDS2666n"))
    w.show()
    w.set_data(data)
    w.handleNewSignals()
    app.exec_()
    w.saveSettings()
