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
#include <QMenu>
#include <QTableView>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <QApplication>
#include <QClipboard>


QPoint translateCoordinate(float x, float y)
{
	return QPoint((int)(1800*x+0.5)+300, 300 - (int)(1800*y+.5));
}

viewport::viewport(int width, int height, const char *title) : scale(width, 70, QImage::Format_RGB32), paintbuff(width, height, QImage::Format_RGB32)
{
	this->setWindowTitle(title);
	this->setFixedSize(width, height+scale.height());
	this->paintbuff.fill(qRgb(255, 255, 255));
	// prepare scale
	this->scale.fill(qRgb(255, 255, 255));
	
	QLinearGradient color(20,0,this->width()-20,0);
	
	color.setColorAt(0, getColorForLevel(mincond));	
	for(float t=0.125f; t < 1.0f; t+= 0.125f) {
	    float level = maxcond*t + (1-t)*mincond;
	    color.setColorAt(t, getColorForLevel(level));	
	}
	color.setColorAt(1, getColorForLevel(maxcond));
	QPainter painter(&scale);
	painter.setBrush(color);
	//painter.setPen(Qt::NoPen);
	painter.drawRect(20, 10, scale.width()-40, 30);	
	painter.setPen(Qt::SolidLine);
	for(int i = 0; i <= 4; i++) {
	  int x = (int)((scale.width()-40)*((float)i/4)+0.5f)+20;
	  float level = (float)i/4;
	  level = maxcond*level + (1-level)*mincond;
	  painter.drawLine(x,40,x,45);
	  painter.drawText(x-25, 45, 51, 20, Qt::AlignHCenter | Qt::AlignTop, QString::number(level));
	  
	}
	  
	
}

QColor viewport::getColorForLevel(float level)
{
	if(level > maxcond) level=maxcond;
	if(level < mincond) level=mincond;
	float fval = ((level-mincond)/(maxcond-mincond));
	fval = std::pow(fval, 1.8f)*255;
	int val = (int)fval;
	return QColor(val, val, val);
}

QBrush viewport::getBrushForElement(float *solution, int n1, int n2, int n3)
{
	// Reorder the elements so that sigma(n1) <= sigma(n2) <= sigma(n3)
	float s1, s2, s3, f_aux;
	int aux;
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
		int level = 255*(((s1+s3)/2));
		return QBrush(getColorForLevel((s1+s3)/2));
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
	
	color.setColorAt(0, getColorForLevel(s1));	
	for(float t=0.125f; t < 1.0f; t+= 0.125f) {
	    float level = s3*t + (1-t)*s1;
	    color.setColorAt(t, getColorForLevel(level));	
	}
	color.setColorAt(1, getColorForLevel(s3));

	return QBrush(color);
}

void viewport::solution_updated(const QModelIndex & topLeft, const QModelIndex & bottomRight )
{
      {		// Redraw on the buffer
	    float *solution = (float *)topLeft.internalPointer() - topLeft.row();
	    int i;
	    QPainter painter(&this->paintbuff);
	    painter.setRenderHint(QPainter::Antialiasing);
	    painter.eraseRect(this->geometry());
	    painter.setPen(Qt::NoPen);
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
		    painter.setBrush(getBrushForElement(solution, elements[i].n1, elements[i].n2, elements[i].n3));
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
	    painter.setPen(Qt::SolidLine);
	    painter.setBrush(Qt::NoBrush);
	    painter.drawEllipse(QRect(translateCoordinate(-0.150f, 0.150f), translateCoordinate(0.150f, -0.150f)));
	}
	
	    
	// Request an update
	this->update();
}

void viewport::paintEvent ( QPaintEvent * event )
{
	QPainter painter(this);
	painter.drawImage(event->rect(), this->paintbuff, event->rect());
	
	painter.drawImage(event->rect(), this->scale, event->rect().translated(0,-this->paintbuff.height()));
	
}

TableViewCopyDataPopupMenu::TableViewCopyDataPopupMenu() {}

TableViewCopyDataPopupMenu *TableViewCopyDataPopupMenu::instance = NULL;

TableViewCopyDataPopupMenu *TableViewCopyDataPopupMenu::getInstance() {
	if(instance==NULL) instance = new TableViewCopyDataPopupMenu();
	return instance;
}

/*void TableViewCopyDataPopupMenu::showMenu(const QPoint & pos )
{
      // Verifies if the sender is actually a QTableView
      if(this->sender()->inherits("QTableView")) {
		QMenu newMenu;
		newMenu.addAction("Copy", this, SLOT(actionFired()))->setData(QVariant::fromValue(sender()));
		newMenu.exec(pos);
      }
}*/

void TableViewCopyDataPopupMenu::actionFired()
{
	    QAction *act = (QAction *) this->sender();
	    QTableView *view = dynamic_cast<QTableView *>(act->parent());
	    QModelIndexList selection = view->selectionModel()->selectedIndexes();
	    QString data;
	    QModelIndex current;
	    foreach(current, selection) {
		data.append(current.data(Qt::DisplayRole).toString());
		data.append(" ");
	    }
	    QApplication::clipboard()->setText(data);
}

solutionView::solutionView(int rows) : rows(rows)
{
	this->sol = new float[rows];
	std::fill(this->sol, this->sol+rows, 2.0);
}

solutionView::~solutionView()
{
	delete[] this->sol;
}

int solutionView::rowCount( const QModelIndex & parent ) const
{
	return this->rows;
}

QVariant  solutionView::data( const QModelIndex & index, int role ) const
{
    if (!index.isValid() || role != Qt::DisplayRole)
        return QVariant();
    return this->sol[index.row()];
}


QVariant  solutionView::headerData (
		int section,
		Qt::Orientation  orientation,
		int role) const
 {
     if (role == Qt::DisplayRole)
         return QString::number(section);
     return QVariant();
 }

void solutionView::setCurrentSolution(const float *newsol)
{
	{
		QMutexLocker(&this->solutionMutex);
		std::copy(newsol, newsol + this->rows, this->sol);
	}
	emit dataChanged(this->createIndex(0,0, this->sol), this->createIndex(0,this->rows));
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
 
 
matrix2ViewModel::matrix2ViewModel(const matrix2 &m) : innermatrix(m)
{}

int matrix2ViewModel::rowCount( const QModelIndex & parent ) const
{
        return this->innermatrix.rows();
}

int matrix2ViewModel::columnCount( const QModelIndex & parent ) const
{
        return this->innermatrix.cols();
}

QVariant  matrix2ViewModel::data( const QModelIndex & index, int role ) const
{
    if (!index.isValid() || role != Qt::DisplayRole)
        return QVariant();
    return this->innermatrix.coeff(index.row(), index.column());
}


QVariant  matrix2ViewModel::headerData (
                int section,
                Qt::Orientation  orientation,
                int role) const
 {
     if (role == Qt::DisplayRole)
         return QString::number(section);
     return QVariant();
 }
 
#include "graphics.moc"

