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
#include <string>


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

class viewport : public QWidget {
 Q_OBJECT
		QImage paintbuff;
		float *solution;
		float *solvar;
		QMutex solutionMutex;
		QColor getColorForIndex(int i);
	public:
		viewport(int width, int height, const char *title);
		QImage &getBuffer() {
			return this->paintbuff;
		}
		void setCurrentSolution(float *val, float *var);
		void saveImage(const std::string &filename);
	protected slots:
		void _saveImage(const QString &filename);
	protected:
	  // override default paint event
	  void paintEvent ( QPaintEvent * event );
};

void drawElements(viewport &view);




#endif /* GRAPHICS_H_ */
