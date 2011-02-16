#include <QtGui/QApplication>
#include "eitannealingtest.h"


int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    eitannealingtest foo;
    foo.show();
    return app.exec();
}
