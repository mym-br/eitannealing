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
#include "solver.h"


class matrixViewModel : public QAbstractTableModel {
	private:
		const matrix &innermatrix;
	public:
		matrixViewModel(const matrix &);
		int rowCount( const QModelIndex & parent = QModelIndex() ) const;
		int columnCount( const QModelIndex & parent = QModelIndex() ) const;
		QVariant  data( const QModelIndex & index, int role = Qt::DisplayRole ) const;
		QVariant  headerData ( int section, Qt::Orientation  orientation, int role = Qt::DisplayRole ) const;
};

class solutionView : public QAbstractListModel {
Q_OBJECT
	private:
		float *sol;
		int rows;
		QMutex solutionMutex;
	public:
		solutionView(int rows);
		~solutionView();
		int rowCount( const QModelIndex & parent = QModelIndex() ) const;
		QVariant  data( const QModelIndex & index, int role = Qt::DisplayRole ) const;
		QVariant  headerData ( int section, Qt::Orientation  orientation, int role = Qt::DisplayRole ) const;
		void setCurrentSolution(const float *newsol);	
};

class viewport : public QWidget {
Q_OBJECT
	QImage paintbuff;
	QImage scale;
	public:
	    viewport(int width, int height, const char *title);
	    QImage &getBuffer() {
		    return this->paintbuff;
	    }
	      
      public slots:
	      void solution_updated(const QModelIndex & topLeft, const QModelIndex & bottomRight );
	protected:
	      QBrush getBrushForElement(float *solution, int n1, int n2, int n3);
	      QColor getColorForLevel(float level);
	      // override default paint event
	      void paintEvent ( QPaintEvent * event );
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
