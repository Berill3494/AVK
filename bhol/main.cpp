#include <ncurses.h>
#include <chrono>
#include <thread>

using namespace std;
using namespace chrono_literals;

struct Vector {
	int x, y;

	void add(Vector const& other) {
		this->x += other.x;
		this->y += other.y;
	}
};

struct Particle {
	Vector	p;
	Vector  v;
	char	c;

	void draw(WINDOW *w){
		wmove(w, this->p.y, this->p.x);
		waddch(w, c);
	}

	void move() {
		p.add(v);	
	}

};


struct Screen {
	WINDOW *w = 0;
	size_t x, y = 0;

	void init() {
		initscr();
		curs_set(0);
		getmaxyx(stdscr, this->y, this->x);
		this->w = newwin(y, x, 0, 0);
	}

	void refresh() {
		wrefresh(this->w);
	}
};



int main() {
    Screen s;
	s.init();

	std::vector<Particle> ps{
		{{0, 0}, {1, 1}, 'B'},
		{{1, 0}, {1, 2}, 'E'},
		{{2, 0}, {2, 1}, 'R'},
		{{3, 0}, {4, 2}, 'I'},
		{{4, 0}, {2, 3}, 'L'}};
	for(;;) {
		wclear(s.w);
		for(auto &p: ps) {
			p.move();
			p.draw(s.w);
		}
		s.refresh();
		std::this_thread::sleep_for(100ms);

		for(auto &p: ps) {
			if(p.p.x >= s.x || p.p.x <= 0) p.v.x*=-1;
			if(p.p.y >= s.y || p.p.y <= 0) p.v.y*=-1;
		}

	}
	endwin();
}
