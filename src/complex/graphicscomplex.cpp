/*
 * graphics.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: thiago
 */

#include "graphics.h"
#include "twodim/problem2D.h"
//#include "problemdescription.h"
#include <QColor>
#include <QPolygon>
#include <QPainter>
#include <QPaintEvent>
#include <QMutexLocker>
#include <QMenu>
#include <QTableView>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <QApplication>
#include <QClipboard>


QPoint viewportcomplex::translateCoordinate(float x, float y)
{
	return QPoint((int)(1800*x+0.5)+300, 300 - (int)(1800*y+.5));
}

viewportcomplex::viewportcomplex(int width, int height, const char *title, std::shared_ptr<problem2D<Complex, Eigen::VectorXcd, matrixcomplex>> _input, double mincond, double maxcond) :
scale(width, 70, QImage::Format_RGB32), paintbuff(width, height, QImage::Format_RGB32), input(_input), minval(mincond), maxval(maxcond)
{
	this->setWindowTitle(title);
	this->setFixedSize(width, height+scale.height());
	this->paintbuff.fill(qRgb(255, 255, 255));
	// prepare scale
	this->scale.fill(qRgb(255, 255, 255));
	
	QLinearGradient color(20,0,this->width()-20,0);
	
	color.setColorAt(0, getColorForLevel(minval));
	for(double t=0.125f; t < 1.0f; t+= 0.125f) {
	    double level = maxval*t + (1-t)*minval;
	    color.setColorAt(t, getColorForLevel(level));	
	}
	color.setColorAt(1, getColorForLevel(maxval));
	QPainter painter(&scale);
	painter.setBrush(color);
	//painter.setPen(Qt::NoPen);
	painter.drawRect(20, 10, scale.width()-40, 30);	
	painter.setPen(Qt::SolidLine);
	for(int i = 0; i <= 4; i++) {
	  int x = (int)((scale.width()-40)*((float)i/4)+0.5f)+20;
	  double level = (float)i/4;
	  level = maxval*level + (1-level)*minval;
	  painter.drawLine(x,40,x,45);
	  painter.drawText(x-25, 45, 51, 20, Qt::AlignHCenter | Qt::AlignTop, QString::number(level));
	  
	}
	  
	
}

QColor viewportcomplex::getColorForLevel(double level)
{
	if(level > maxval) level=maxval;
	if(level < minval) level=minval;
	double fval = ((level-minval)/(maxval-minval));
	fval = std::pow(fval, 1.8f)*255;
	int val = (int)fval;
	return QColor(val, val, val);
}

QBrush viewportcomplex::getBrushForElement(double *solution, int n1, int n2, int n3)
{
	// Reorder the elements so that sigma(n1) <= sigma(n2) <= sigma(n3)
	double s1, s2, s3, f_aux;
	int aux;
	s1 = solution[input->node2coefficient[n1]];
	s2 = solution[input->node2coefficient[n2]];
	s3 = solution[input->node2coefficient[n3]];
	// Quick hardcoded sort
	
	if(s1 > s2) {
	    f_aux = s1; s1 = s2; s2 = f_aux;
	    aux = n1; n1 = n2; n2 = aux;
	}	
	if(s2 > s3) {
	    f_aux = s2; s2 = s3; s3 = f_aux;
	    aux = n2; n2 = n3; n3 = aux;
	    if(s1 > s2) {
		f_aux = s1; s1 = s2; s2 = f_aux;
		aux = n1; n1 = n2; n2 = aux;
	    }
	}	
	if((s3-s1)<0.001) {
		int level = 255*(((s1+s3)/2));
		return QBrush(getColorForLevel((s1+s3)/2));
	}
	
	// Ok, so now we must find both control points
	// 	let's say p0 = n1		
	double alpha = (s2-s1)/(s3-s1);	// 0 <= x <= 1
	Eigen::Vector2d a(
		input->nodes[n3].x - input->nodes[n1].x,
		input->nodes[n3].y - input->nodes[n1].y),
					b(
		input->nodes[n2].x - input->nodes[n1].x,
		input->nodes[n2].y - input->nodes[n1].y);

	double c = b.y() - alpha*a.y();
	double s = alpha*a.x() - b.x();
	double cc = c*c;
	double ss = s*s;
	double cs = c*s;
	double x = cc*a.x() + cs*a.y();
	double y = ss*a.y() + cs*a.x();
	x /= (cc+ss);
	y /= (cc+ss);

	QLinearGradient color(
		translateCoordinate(input->nodes[n1].x, input->nodes[n1].y),
		translateCoordinate(input->nodes[n1].x + x, input->nodes[n1].y + y));
	
	color.setColorAt(0, getColorForLevel(s1));	
	for(double t=0.125f; t < 1.0f; t+= 0.125f) {
	    double level = s3*t + (1-t)*s1;
	    color.setColorAt(t, getColorForLevel(level));	
	}
	color.setColorAt(1, getColorForLevel(s3));

	return QBrush(color);
}

void viewportcomplex::solution_updated(const QModelIndex & topLeft, const QModelIndex & bottomRight)
{
      {		// Redraw on the buffer
	    double *solution = (double *)topLeft.internalPointer() - topLeft.row();
	    int i;
	    QPainter painter(&this->paintbuff);
	    painter.setRenderHint(QPainter::Antialiasing);
	    painter.eraseRect(this->geometry());
	    painter.setPen(Qt::NoPen);
		for (i = 0; i<input->elements.size(); i++) {
		    QPolygon polygon;
			int n = input->elements[i].a;
			polygon.push_back(translateCoordinate(input->nodes[n].x, input->nodes[n].y));
			n = input->elements[i].b;
			polygon.push_back(translateCoordinate(input->nodes[n].x, input->nodes[n].y));
			n = input->elements[i].c;
			polygon.push_back(translateCoordinate(input->nodes[n].x, input->nodes[n].y));
			n = input->elements[i].a;
			polygon.push_back(translateCoordinate(input->nodes[n].x, input->nodes[n].y));
		    int level = 255;//*(solution[elements[i].condIndex]-1);		
			painter.setBrush(getBrushForElement(solution, input->elements[i].a, input->elements[i].b, input->elements[i].c));
		    painter.drawConvexPolygon(polygon);
	    }
	    // Draw electrodes
	    /*for(i=0;i<electrodes.size();i++) {
		    node *base = &nodes[electrodes[i].baseNode];
		    node *n1 = &nodes[electrodes[i].n1];
		    node *n2 = &nodes[electrodes[i].n2];
		    node *n3 = &nodes[electrodes[i].n3];
		    painter.drawLine(
				    translateCoordinate(base->x, base->y),
				    translateCoordinate(n1->x, n1->y));
		    painter.drawLine(
				    translateCoordinate(base->x, base->y),
				    translateCoordinate(n2->x, n2->y));
		    painter.drawLine(
				    translateCoordinate(base->x, base->y),
				    translateCoordinate(n3->x, n3->y));
	    }*/
	    
	    // Draw perimeter
	    
	    painter.setPen(Qt::SolidLine);
	    painter.setBrush(Qt::NoBrush);
		for (auto const &edge : input->perimeter) {
	      painter.drawLine(
			  translateCoordinate(input->nodes[edge.first].x, input->nodes[edge.first].y),
			  translateCoordinate(input->nodes[edge.second].x, input->nodes[edge.second].y)
	      );
	      
	    }
	}
	
	    
	// Request an update
	this->update();
}

void viewportcomplex::paintEvent(QPaintEvent * event)
{
	QPainter painter(this);
	painter.drawImage(event->rect(), this->paintbuff, event->rect());
	
	painter.drawImage(event->rect(), this->scale, event->rect().translated(0,-this->paintbuff.height()));
	
}