/*
 * graphics.h
 *
 *  Created on: Jun 25, 2010
 *      Author: thiago
 */

#ifndef GRAPHICS_H_
#define GRAPHICS_H_

#include <QApplication>
#include <QImage>
#include <QWidget>
#include <QPainter>
#include <QAbstractTableModel>
#include <QMutex>
#include <memory>
#include "twodim\problem2D.h"
//class problem2D;

class solutionView : public QAbstractListModel {
Q_OBJECT
	private:
		double *sol;
		int rows;
		QMutex solutionMutex;
	public:
		solutionView(int rows);
		~solutionView();
		int rowCount( const QModelIndex & parent = QModelIndex() ) const;
		QVariant  data( const QModelIndex & index, int role = Qt::DisplayRole ) const;
		QVariant  headerData ( int section, Qt::Orientation  orientation, int role = Qt::DisplayRole ) const;
		void setCurrentSolution(const double *newsol);	
};

class viewport : public QWidget {
Q_OBJECT
	public:
		viewport(int width, int height, const char *title, std::shared_ptr<problem2D> _input);
	    QImage &getBuffer() {
		    return this->paintbuff;
	    }
	      
      public slots:
		void solution_updated(const QModelIndex & topLeft, const QModelIndex & bottomRight);
	protected:
		  QImage paintbuff;
		  QImage scale;
	      QBrush getBrushForElement(double *solution, int n1, int n2, int n3);
	      QColor getColorForLevel(double level);
	      // override default paint event
	      void paintEvent ( QPaintEvent * event );
		  std::shared_ptr<problem2D> input;
};

class viewportcomplex : public QWidget {
	Q_OBJECT
public:
	viewportcomplex(int width, int height, const char *title, std::shared_ptr<problem2D> _input, double mincond, double maxcond);
	QImage &getBuffer() {
		return this->paintbuff;
	}

	public slots:
	void solution_updated(const QModelIndex & topLeft, const QModelIndex & bottomRight);
protected:
	QImage paintbuff;
	QImage scale;
	QBrush getBrushForElement(double *solution, int n1, int n2, int n3);
	QColor getColorForLevel(double level);
	// override default paint event
	void paintEvent(QPaintEvent * event);
	std::shared_ptr<problem2D> input;
private:
	QPoint translateCoordinate(float x, float y);
	double minval, maxval;
};

class TableViewCopyDataPopupMenu : public QObject{
Q_OBJECT
      private:
	    // Singleton
	    static TableViewCopyDataPopupMenu *instance;
	    TableViewCopyDataPopupMenu();
      public:
	    static TableViewCopyDataPopupMenu *getInstance();
      public slots:
	    void actionFired();  
};

#endif /* GRAPHICS_H_ */
