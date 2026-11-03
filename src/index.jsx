import React from 'react';
import ReactDOM from 'react-dom/client';
import './styles/explorer.css';
import App from './App';

// Theme before the first paint: the saved choice, else the system's
let theme = null;
try { theme = localStorage.getItem('ascatsc-theme'); } catch (e) { /* private mode */ }
document.documentElement.dataset.theme = theme || (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');

const root = ReactDOM.createRoot(document.getElementById('root'));
root.render(
  <React.StrictMode>
    <App />
  </React.StrictMode>
);
