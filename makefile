all: package

package:
	zip popcodeSuite.zip *.R *.html *.css

deploy:
	mkdir -p $${HOME}/www/html/popcodeSuite
	mkdir -p $${HOME}/www/cgi/popcodeSuite
	chmod o+w $${HOME}/www/html/popcodeSuite
	ln -s cgiLauncher.R $${HOME}/www/cgi/popcodeSuite/cgiLauncher.R
	ln -s cgiInput.html $${HOME}/www/html/popcodeSuite/cgiInput.html
	ln -s error.html $${HOME}/www/html/popcodeSuite/error.html
	ln -s result.html $${HOME}/www/html/popcodeSuite/cgiLauncher.html
	ln -s style.css $${HOME}/www/html/popcodeSuite/style.css
	ln -s wait.html $${HOME}/www/html/popcodeSuite/wait.html