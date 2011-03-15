/*
 * graphics.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: thiago
 */

#include "graphics.h"
#include "problemdescription.h"
#include <QColor>
#include <QPolygon>
#include <QPainter>
#include <QPaintEvent>
#include <QMutexLocker>
#include <Eigen/Core>
#include <Eigen/Geometry>


QPoint translateCoordinate(float x, float y)
{
	return QPoint((int)(1800*x+0.5)+300, 300 - (int)(1800*y+.5));
}

viewport::viewport(int width, int height, const char *title) : paintbuff(width, height, QImage::Format_RGB32)
{
	this->solution = new float[numcoefficients];
	for(int i=0;i<numcoefficients;i++) solution[i] = 2;
	this->setWindowTitle(title);
	this->setFixedSize(width, height);
	QPainter(this).fillRect(0, 0, width, height, qRgb(255, 255, 255));
	this->paintbuff.fill(qRgb(255, 255, 255));
}

QBrush viewport::getBrushForElement(int n1, int n2, int n3)
{
	// Reorder the elements so that sigma(n1) <= sigma(n2) <= sigma(n3)
	float s1, s2, s3;
	int aux, f_aux;
	s1 = solution[node2coefficient[n1]];
	s2 = solution[node2coefficient[n2]];
	s3 = solution[node2coefficient[n3]];
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
		int level = 255*(((s1+s3)/2)-1);
		return QBrush(QColor(level,level,level));
	}
	
	// Ok, so now we must find both control points
	// 	let's say p0 = n1		
	float alpha = (s2-s1)/(s3-s1);	// 0 <= x <= 1
	Eigen::Vector2d a(
		nodes[n3].x - nodes[n1].x,
		nodes[n3].y - nodes[n1].y),
					b(
		nodes[n2].x - nodes[n1].x,
		nodes[n2].y - nodes[n1].y);

	float c = b.y() - alpha*a.y();
	float s = alpha*a.x() - b.x();
	float cc = c*c;
	float ss = s*s;
	float cs = c*s;
	float x = cc*a.x() + cs*a.y();
	float y = ss*a.y() + cs*a.x();
	x /= (cc+ss);
	y /= (cc+ss);

	QLinearGradient color(
		translateCoordinate(nodes[n1].x, nodes[n1].y),
		translateCoordinate(nodes[n1].x+x, nodes[n1].y+y));
	
	int level = 255*(s1-1);
	color.setColorAt(0, QColor(level,level,level));
	level = 255*(s3-1);
	color.setColorAt(1, QColor(level,level,level));

	return QBrush(color);
}

void drawElements(viewport &view)
{
	int i;
	QPainter painter(&(view.getBuffer()));
	painter.eraseRect(view.geometry());
	for(i=0;i<elements.size();i++) {
		QPolygon polygon;
		int n = elements[i].n1;
		polygon.push_back(translateCoordinate(nodes[n].x, nodes[n].y));
		n = elements[i].n2;
		polygon.push_back(translateCoordinate(nodes[n].x, nodes[n].y));
		n = elements[i].n3;
		polygon.push_back(translateCoordinate(nodes[n].x, nodes[n].y));
		n = elements[i].n1;
		polygon.push_back(translateCoordinate(nodes[n].x, nodes[n].y));
		painter.drawConvexPolygon(polygon);
	}
	// Draw electrodes
	for(i=0;i<electrodes.size();i++) {
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
	}
}

void drawSolution(viewport &view, float *solution)
{

}

void viewport::setCurrentSolution(float *val)
{
	QMutexLocker lock(&solutionMutex);
	memcpy(this->solution, val, sizeof(float)*numcoefficients);
	QMetaObject::invokeMethod(this, "update", Qt::QueuedConnection);
}

void viewport::paintEvent ( QPaintEvent * event )
{
	//QPainter(this).drawImage(event->rect(), this->paintbuff, event->rect());
	QMutexLocker lock(&solutionMutex);
	int i;
	QPainter painter(this);
	painter.setRenderHint(QPainter::Antialiasing);
	painter.eraseRect(this->geometry());
	for(i=0;i<elements.size();i++) {
		QPolygon polygon;
		int n = elements[i].n1;
		polygon.push_back(translateCoordinate(nodes[n].x, nodes[n].y));
		n = elements[i].n2;
		polygon.push_back(translateCoordinate(nodes[n].x, nodes[n].y));
		n = elements[i].n3;
		polygon.push_back(translateCoordinate(nodes[n].x, nodes[n].y));
		n = elements[i].n1;
		polygon.push_back(translateCoordinate(nodes[n].x, nodes[n].y));
		int level = 255;//*(solution[elements[i].condIndex]-1);		
		painter.setBrush(getBrushForElement(elements[i].n1, elements[i].n2, elements[i].n3));
		painter.drawConvexPolygon(polygon);
	}
	// Draw electrodes
	for(i=0;i<electrodes.size();i++) {
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
	}
}

matrixViewModel::matrixViewModel(const matrix &m) : innermatrix(m)
{}

int matrixViewModel::rowCount( const QModelIndex & parent ) const
{
	return this->innermatrix.rows();
}

int matrixViewModel::columnCount( const QModelIndex & parent ) const
{
	return this->innermatrix.cols();
}

QVariant  matrixViewModel::data( const QModelIndex & index, int role ) const
{
    if (!index.isValid() || role != Qt::DisplayRole)
        return QVariant();
    if(index.column()<index.row())
    	return this->innermatrix.coeff(index.row(), index.column());
	return this->innermatrix.coeff(index.column(), index.row());
}


QVariant  matrixViewModel::headerData (
		int section,
		Qt::Orientation  orientation,
		int role) const
 {
     if (role == Qt::DisplayRole)
         return QString::number(section);
     return QVariant();
 }


